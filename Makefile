CODE_DIR  = ~/phantom-grtest
BUILD_DIR = $(CODE_DIR)/build
CODE_NAME = grtest

default:
	cd $(BUILD_DIR) && make
	cp $(BUILD_DIR)/$(CODE_NAME) .

clean:
	cd $(BUILD_DIR) && make clean
