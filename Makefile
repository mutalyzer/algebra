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

CC         = gcc
CFLAGS     = -std=c11 -march=native -Wall -Wextra -Wpedantic -fanalyzer

MEMCHECK   = valgrind -q --leak-check=full --error-exitcode=1


.PHONY: check clean debug release


debug:   CFLAGS+=-O0 -g -DDEBUG -mno-avx512f
release: CFLAGS+=-O3


debug release: $(TARGET)


check: CFLAGS+=-O0 -g -DDEBUG -ftest-coverage -fprofile-arcs -dumpbase '' \
		-fkeep-inline-functions -fkeep-static-functions -mno-avx512f
check: $(TEST_BINS) $(TESTS_GCNOS) $(TEST_GCDAS) $(TEST_COVS)


clean:
	$(RM) $(TARGET) $(OBJECTS) $(DEPS) \
		$(TEST_BINS) $(TEST_DEPS) $(TEST_GCNOS) $(TEST_GCDAS) \
		$(TEST_DIR)/*.gcov


$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^


-include $(TEST_DEPS)

%.out: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -o $@ $<

%.gcda: %.out
	$(MEMCHECK) ./$<

%.c.gcov: %.gcda
	gcov $<
	mv *.gcov $(TEST_DIR)


-include $(DEPS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -o $@ -c $<
