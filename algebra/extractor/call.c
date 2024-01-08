#include <stddef.h>     // size_t
#include <stdio.h>      // printf
#include <stdlib.h>     // EXIT_*
#include <string.h>     // str*


struct PMR
{
    size_t start;
    size_t period;
    size_t count;
    size_t remainder;
}; // PMR


size_t find_runs(char const* word, struct PMR* pmrs);


static inline size_t
max(size_t const lhs, size_t const rhs)
{
    return lhs > rhs ? lhs : rhs;
} // max


static char const*
nth_fibonacci_word(size_t const n) // als n < 2 dragons
{
    size_t last_len = 1;
    char* last_word = malloc(last_len + 1);
    last_word[0] = 'C';
    last_word[1] = '\0';

    size_t current_len = 1;
    char* current_word = malloc(current_len + 1);
    current_word[0] = 'A';
    current_word[1] = '\0';

    for (size_t i = 2; i < n; i++)
    {
        size_t const len = last_len + current_len;
        char* next_word = malloc(len + 1);
        if (next_word == NULL)
        {
            printf("malloc()\n");
            return NULL;
        } // if
        strncpy(next_word, current_word, current_len);
        strncpy(next_word + current_len, last_word, last_len);
        next_word[len] = '\0';

        free(last_word);

        last_word = current_word;
        last_len = current_len;
        current_word = next_word;
        current_len = len;
    } // for

    free(last_word);
    return current_word;
} // fibo


int
main(void)
{
    char const* word = nth_fibonacci_word(42);
    size_t const len = strlen(word);

    struct PMR* pmrs = malloc(sizeof(*pmrs) * len);
    if (pmrs == NULL)
    {
        printf("malloc()\n");
        return EXIT_FAILURE;
    } // if

    size_t const n = find_runs(word, pmrs);

    size_t* depth = malloc(sizeof(*depth) * len);
    for (size_t i = 0; i < len; ++i)
    {
        depth[i] = 0;
    } // for

    size_t sum = 0;
    size_t max_depth = 0;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t p = pmrs[i].period * 2 - 1; p < pmrs[i].period * pmrs[i].count + pmrs[i].remainder; ++p)
        {
            depth[pmrs[i].start + p] += 1;
            max_depth = max(max_depth, depth[pmrs[i].start + p]);
            sum += 1;
        } // for
    } // for

    printf("%zu %zu %zu %zu\n", len, n, max_depth, sum);

    free(pmrs);
    free((char*) word);

    return EXIT_SUCCESS;
} // main

