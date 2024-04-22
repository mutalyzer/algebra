SRC_DIR    = src
TEST_DIR   = tests

SOURCES    = $(wildcard $(SRC_DIR)/*.c)
OBJECTS    = $(SOURCES:.c=.o)
DEPS       = $(SOURCES:.c=.d)
TESTS      = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS  = $(TESTS:.c=.out)
TEST_DEPS  = $(TESTS:.c=.d)
TEST_GCNOS = $(TESTS:.c=.gcno)
TEST_GCDAS = $(TESTS:.c=.gcda)
TEST_COVS  = $(TESTS:.c=.c.gcov)

TARGET     = a.out

CC         = gcc-11
CFLAGS     = -std=c99 -march=native -Wall -Wextra -Wpedantic -fanalyzer

MEMCHECK   = valgrind -q --leak-check=full --error-exitcode=1


.PHONY: check clean debug release


debug:   CFLAGS+=-O0 -g -DDEBUG
release: CFLAGS+=-O3


debug release: $(TARGET)


check: CFLAGS+=-O0 -g -DDEBUG -ftest-coverage -fprofile-arcs -dumpbase '' \
		-fkeep-inline-functions -fkeep-static-functions
check: $(TEST_BINS) $(TESTS_GCNOS) $(TEST_GCDAS) $(TEST_COVS)


clean:
	$(RM) $(TARGET) $(OBJECTS) $(DEPS) \
		$(TEST_BINS) $(TEST_DEPS) $(TEST_GCNOS) $(TEST_GCDAS) \
		$(TEST_DIR)/*.c.gcov


$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^


-include $(TEST_DEPS)

%.out: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -o $@ $<

%.gcda: %.out
	$(MEMCHECK) ./$<

%.c.gcov: %.gcda
	gcov-11 $<
	mv *.c.gcov $(TEST_DIR)


-include $(DEPS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -o $@ -c $<
