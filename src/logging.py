import sys
import logging
import time
import os

logger = logging.getLogger("python_ibd")


class StageLogger:
    """Logs a stage with an inline progress bar on stdout and DEBUG messages to the log file.

    All stdout prints and log file entries are optionally prefixed with stage_name.
    Call StageLogger.setup(output_prefix) once before creating any instances to configure the log file.

    Usage:
        StageLogger.setup("output/my_run")  # creates output/my_run.log
        log = StageLogger("Stage 2/5")
        log.stdout("Computing Allele Frequencies...")
        log.debug("Started computing allele frequencies for 2000 variants...")
        for i in range(total):
            # ... work ...
            log.update_progress(i, total)
        log.finish("Done.")
    """

    @staticmethod
    def setup(output_prefix: str):
        """Configure the log file at {output_prefix}.log. Call once before creating StageLogger instances."""
        filepath = f"{output_prefix}.log"
        dir_name = os.path.dirname(filepath)
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)

        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        file_fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        fh = logging.FileHandler(filepath, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(file_fmt)
        logger.addHandler(fh)

    def __init__(
        self,
        stage_name: str = None,
        bar_width: int = 30,
        min_redraw_interval: float = 0.033,
    ):
        self._stage_name = stage_name
        self._prefix = f"{stage_name} " if stage_name else ""
        self._log_prefix = f"{stage_name} - " if stage_name else ""
        self._bar_width = bar_width
        self._min_redraw_interval = min_redraw_interval
        self._last_draw_time = 0.0

    def stdout(self, text: str):
        logger.debug("%s%s", self._log_prefix, text)
        sys.stdout.write(f"\r{text}\033[K")
        sys.stdout.flush()
        self.last_stdout_text = text

    def debug(self, text: str):
        logger.debug("%s%s", self._log_prefix, text)

    def info(self, text: str):
        logger.info("%s%s", self._log_prefix, text)

    def warning(self, text: str):
        logger.warning("%s%s", self._log_prefix, text)

    def error(self, text: str):
        logger.error("%s%s", self._log_prefix, text)

    def update_progress(self, current: int, total: int, message: str = "Progress"):
        is_last = current >= total - 1
        now = time.time()
        if not is_last and (now - self._last_draw_time) < self._min_redraw_interval:
            return
        self._last_draw_time = now

        pct = (current + 1) / total
        filled = int(self._bar_width * pct)
        bar = (
            "=" * filled
            + (">" if filled < self._bar_width else "")
            + " " * max(0, self._bar_width - filled - 1)
        )
        sys.stdout.write(
            f"\r{self.last_stdout_text} | {message} [{bar}] {pct * 100:5.1f}%\033[K"
        )
        sys.stdout.flush()

    def finish(self, text: str):
        logger.debug("%s%s", self._log_prefix, text)
        sys.stdout.write(f"\r{text}\033[K\n")
        sys.stdout.flush()
