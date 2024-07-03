# Location of top-level MicroPython directory
MPY_DIR = ../../..

# Name of module
MOD = thumbyrt

# Source files (.c or .py)
SRC = main.c qfplib.S lmul.S ldivmod.S uldivmod.S

# Architecture to build for (x86, x64, armv7m, xtensa, xtensawin)
ARCH = armv6m

# Include to get the rules for compiling and linking the module
include $(MPY_DIR)/py/dynruntime.mk
