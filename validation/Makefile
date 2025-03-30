CC = gcc
CFLAGS = -Wall -g -I.
LDFLAGS = -lm

SOURCES = validation.c perf.c problemDef.c util.c
BUILD_DIR = build
SRC_DIR = .
DATA_DIR = $(BUILD_DIR)/data
OBJECTS = $(addprefix $(BUILD_DIR)/, $(notdir $(SOURCES:.c=.o)))
EXECUTABLE = $(BUILD_DIR)/validate

all: $(BUILD_DIR) $(EXECUTABLE) $(DATA_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(DATA_DIR):
	mkdir -p $(DATA_DIR)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean 