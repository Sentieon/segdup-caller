import colorlog
import logging

# Global log level that can be set programmatically
_global_log_level = None

handler = colorlog.StreamHandler()
handler.setFormatter(
    colorlog.ColoredFormatter("%(log_color)s%(levelname)s:%(name)s:%(message)s")
)
# Set a default level for the handler
handler.setLevel(logging.INFO)


def set_global_log_level(level):
    """Set the global log level for all loggers."""
    global _global_log_level
    if isinstance(level, str):
        _global_log_level = getattr(logging, level.upper())
    else:
        _global_log_level = level
    
    # Also set the handler level to ensure messages are passed through
    handler.setLevel(_global_log_level)
    
    # Update all existing loggers that were created with the old level
    updated_count = 0
    for name in logging.Logger.manager.loggerDict:
        logger_obj = logging.Logger.manager.loggerDict[name]
        
        # Some entries in loggerDict might be PlaceHolder objects, not actual loggers
        if isinstance(logger_obj, logging.Logger):
            if hasattr(logger_obj, 'handlers'):
                # Check if any of the logger's handlers is our shared handler
                for logger_handler in logger_obj.handlers:
                    if logger_handler is handler:
                        logger_obj.setLevel(_global_log_level)
                        updated_count += 1
                        break
    
    # Also update the handler level for any logger that uses our shared handler
    # This is a more direct approach
    if hasattr(handler, 'loggers_using_this_handler'):
        for logger_ref in handler.loggers_using_this_handler:
            if logger_ref() is not None:  # WeakRef check
                logger_ref().setLevel(_global_log_level)


def get_logger(name: str, level=None):
    """Return a logger with a colorlog handler."""
    logger = colorlog.getLogger(name)
    
    # Only add handler if not already added (to avoid duplicates)
    if handler not in logger.handlers:
        logger.addHandler(handler)
    
    logger.propagate = False
    
    # Use global log level if set and no specific level provided, default to INFO
    effective_level = level or _global_log_level or logging.INFO
    if isinstance(effective_level, str):
        effective_level = getattr(logging, effective_level.upper())
    logger.setLevel(effective_level)
    
    return logger
