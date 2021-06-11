from typing import List, Optional
from pytorch_lightning import LightningModule, LightningDataModule, Callback, Trainer
from pytorch_lightning.loggers import LightningLoggerBase
from pytorch_lightning import seed_everything
import hydra
from omegaconf import DictConfig
from src.utils import utils
from src.utils.cv.data_module import KFoldCVDataModule


log = utils.get_logger(__name__)


def train_cv(config: DictConfig) -> Optional[float]:
    """Contains training pipeline.
    Instantiates all PyTorch Lightning objects from config.

    Args:
        config (DictConfig): Configuration composed by Hydra.

    Returns:
        Optional[float]: Metric score for hyperparameter optimization.
    """

    # Set seed for random number generators in pytorch, numpy and python.random
    if "seed" in config:
        seed_everything(config.seed)

    n_splits = config.n_splits

    # Init Lightning datamodule
    log.info(f"Instantiating datamodule <{config.datamodule._target_}>")
    datamodule: LightningDataModule = hydra.utils.instantiate(config.datamodule)

    cv_datamodule = KFoldCVDataModule(datamodule, n_splits)

    ckpt_name = config.callbacks.model_checkpoint["filename"]

    for fold_idx, loaders in enumerate(cv_datamodule.split()):

        config.logger.csv["version"] = f"fold_{fold_idx}"
        config.callbacks.model_checkpoint["filename"] = ckpt_name + f"_fold_{fold_idx}"

        config.logger.wandb["version"] = "_${now:%Y-%m-%d_%H-%M-%S}" + f"fold_{fold_idx}"
        config.logger.wandb["project"] = config.project_name

        # Init Lightning model
        log.info(f"Instantiating model <{config.model._target_}>")
        model: LightningModule = hydra.utils.instantiate(config.model)

        # Init Lightning callbacks
        callbacks: List[Callback] = []
        if "callbacks" in config:
            for _, cb_conf in config["callbacks"].items():
                if "_target_" in cb_conf:
                    log.info(f"Instantiating callback <{cb_conf._target_}>")
                    callbacks.append(hydra.utils.instantiate(cb_conf))

        # Init Lightning loggers
        logger: List[LightningLoggerBase] = []
        if "logger" in config:
            for _, lg_conf in config["logger"].items():
                if "_target_" in lg_conf:
                    log.info(f"Instantiating logger <{lg_conf._target_}>")
                    logger.append(hydra.utils.instantiate(lg_conf))

        # Init Lightning trainer
        log.info(f"Instantiating trainer <{config.trainer._target_}>")
        trainer: Trainer = hydra.utils.instantiate(
            config.trainer, callbacks=callbacks, logger=logger, _convert_="partial"
        )

        # Send some parameters from config to all lightning loggers
        log.info("Logging hyperparameters!")
        utils.log_hyperparameters(
            config=config,
            model=model,
            datamodule=datamodule,
            trainer=trainer,
            callbacks=callbacks,
            logger=logger,
        )

        # Train the model
        log.info("Starting training!")
        trainer.fit(model, *loaders)

        # Evaluate model on test set after training
        if not config.trainer.get("fast_dev_run"):
            log.info("Starting testing!")
            test_dataloader = cv_datamodule.get_test_dataloader()
            if len(test_dataloader) > 0:
                trainer.test(model, test_dataloader)
            else:
                log.info("Test data is empty!")

        # Make sure everything closed properly
        log.info("Finalizing!")
        utils.finish(
            config=config,
            model=model,
            datamodule=datamodule,
            trainer=trainer,
            callbacks=callbacks,
            logger=logger,
        )

        # Print path to best checkpoint
        log.info(f"Best checkpoint path:\n{trainer.checkpoint_callback.best_model_path}")

        # Return metric score for hyperparameter optimization
        optimized_metric = config.get("optimized_metric")
        if optimized_metric:
            return trainer.callback_metrics[optimized_metric]
