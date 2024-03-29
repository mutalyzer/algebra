SRC_DIR  = src
TEST_DIR = tests

SOURCES  = $(sort $(shell find $(SRC_DIR) -name '*.c'))
OBJECTS  = $(SOURCES:.c=.o)
DEPS     = $(OBJECTS:.o=.d)

TARGET   = a.out

TESTS    = $(sort $(shell find $(TEST_DIR) -name 'test_*.c'))

CC       = gcc
CFLAGS   = -std=c11 -march=native -Wall -Wextra -Wpedantic -O0 -g -DDEBUG 

.PHONY: all check clean

all: $(TARGET)


check: CFLAGS += --coverage -fprofile-arcs -fprofile-dir=$(TEST_DIR)
check: $(TESTS:.c=.out)
	valgrind -q ./$^
	gcov $(^:.out=.gcda)
	mv *.gcov $(TEST_DIR)


clean:
	$(RM) $(TARGET) $(OBJECTS) $(DEPS) $(TESTS:.c=.out) $(TEST_DIR)/*.gcno $(TEST_DIR)/*.gcda $(TEST_DIR)/*.gcov


$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^

-include $(DEPS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -o $@ -c $<

%.out: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $<
	mv *.gcno $(TEST_DIR)
