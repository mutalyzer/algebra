#include <limits.h>     // CHAR_BIT
#include <stdbool.h>    // bool
#include <stddef.h>     // size_t
#include <stdio.h>      // FILE, stderr, stdin, stdout, fflush, fprintf, fscanf, printf
#include <stdlib.h>     // EXIT_*, atoll, exit, qsort
#include <string.h>     // strlen, strncmp


#define WORD_SIZE (sizeof(size_t) * CHAR_BIT + 1)


static inline size_t
max(size_t const lhs, size_t const rhs)
{
    return lhs > rhs ? lhs : rhs;
} // max


static inline size_t
min(size_t const lhs, size_t const rhs)
{
    return lhs < rhs ? lhs : rhs;
} // min


struct PMR
{
    size_t start;
    size_t period;
    size_t count;
    size_t remainder;
}; // PMR


static int
pmr_cmp(void const* lhs, void const* rhs)
{
    struct PMR const* lhs_pmr = lhs;
    struct PMR const* rhs_pmr = rhs;
    return lhs_pmr->start - rhs_pmr->start;
} // pmr_cmp


static size_t
decode(size_t const len, char word[len], size_t idx)
{
    size_t const n = len - __builtin_clzll(idx + 1) - 2;
    size_t i = n;
    word[0] = 'A';
    while (idx > 0)
    {
        idx -= 1;
        word[i] = "AC"[idx % 2];
        i -= 1;
        idx /= 2;
    } // while
    word[n + 1] = '\0';
    return n + 1;
} // decode


static bool
are_isomorphic_words(size_t const len, char const lhs[len], char const rhs[len])
{
    char lhs_map[256] = {'\0'};
    char rhs_map[256] = {'\0'};

    for (size_t i  = 0; i < len; ++i)
    {
        if (lhs_map[(unsigned) lhs[i]] != rhs_map[(unsigned) rhs[i]])
        {
            return false;
        } // if
        lhs_map[(unsigned) lhs[i]] = rhs[i];
        rhs_map[(unsigned) rhs[i]] = rhs[i];
    } // for
    return true;
} // are_isomorphic_words


static size_t
find_pmrs_alt(size_t const len, char const word[len], struct PMR pmrs[len])
{
    size_t i = 0;
    for (size_t start = 0; start < len; ++start)
    {
        for (size_t period = 1; period <= (len - start) / 2; ++period)
        {
            size_t count = 1;
            while (strncmp(word + start, word + start + count * period, period) == 0)
            {
                count += 1;
            } // while
            if (count < 2)
            {
                continue;
            } // if

            size_t remainder = 0;
            while (word[start + remainder] == word[start + period * count + remainder])
            {
                remainder += 1;
            } // while

            bool uniq = true;
            for (size_t j = 0; j < i; ++j)
            {
                size_t idx = i - j - 1;
                if (period % pmrs[idx].period == 0 &&
                    pmrs[idx].start + pmrs[idx].count * pmrs[idx].period + pmrs[idx].remainder == start + count * period + remainder)
                {
                    uniq = false;
                    break;
                } // if
            } // for

            if (uniq)
            {
                pmrs[i].start = start;
                pmrs[i].period = period;
                pmrs[i].count = count;
                pmrs[i].remainder = remainder;
                i += 1;
            } // if
        } // for
    } // for
    return i;
} // find_pmrs_alt


static size_t
find_pmrs(size_t const len, char const word[len], struct PMR pmrs[len])
{
    size_t inv[len];
    for (size_t i = 0; i < len; ++i)
    {
        inv[i] = -1;
    } // for

    size_t i = 0;
    for (size_t period = 1; period < len / 2 + 1; ++period)
    {
        size_t start = 0;
        while (start < len - period)
        {
            size_t count = 1;
            while (strncmp(word + start, word + start + period * count, period) == 0)
            {
                count += 1;
            } // while

            size_t remainder = 0;
            while (word[start + remainder] == word[start + period * count + remainder])
            {
                remainder += 1;
            } // while

            if (inv[start] != (size_t) -1)
            {
                struct PMR const pmr = pmrs[inv[start]];
                if (pmr.start + pmr.period * pmr.count + pmr.remainder == start + period * count + remainder &&
                    period % pmr.period == 0)
                {
                    start += pmr.period * pmr.count - period + 1;
                    continue;
                } // if
            } // if

            if (count > 1)
            {
                inv[start] = i;
                pmrs[i].start = start;
                pmrs[i].period = period;
                pmrs[i].count = count;
                pmrs[i].remainder = remainder;
                i += 1;
                start += period * (count - 1) + 1;
            } // if
            else
            {
                start += remainder + 1;
            } // else
        } // while
    } // for
    // qsort(pmrs, i, sizeof(*pmrs), pmr_cmp);
    return i;
} // find_pmrs


static size_t
inv_array(size_t const len, size_t const n, struct PMR const pmrs[n], size_t inv[len][n + 1])
{
    for (size_t i = 0; i < len; ++i)
    {
        inv[i][0] = 0;
    } // for

    size_t depth = 0;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t p = pmrs[i].period * 2 - 1; p < pmrs[i].period * pmrs[i].count + pmrs[i].remainder; ++p)
        {
            size_t const len = inv[pmrs[i].start + p][0];
            inv[pmrs[i].start + p][len + 1] = i;
            inv[pmrs[i].start + p][0] += 1;
            depth = max(depth, len + 1);
        } // for
    } // for
    return depth;
} // inv_array


static void
overlap_array(size_t const n, struct PMR const pmrs[n], size_t overlap[n])
{
    overlap[0] = 0;
    for (size_t rhs = 1; rhs < n; ++rhs)
    {
        overlap[rhs] = 0;
        for (size_t lhs = 0; lhs < rhs; ++lhs)
        {
            size_t const lhs_end = pmrs[lhs].start + pmrs[lhs].period * pmrs[lhs].count + pmrs[lhs].remainder;
            size_t const rhs_end = pmrs[rhs].start + pmrs[rhs].period * pmrs[rhs].count + pmrs[rhs].remainder;
            if (pmrs[lhs].start < pmrs[rhs].start && pmrs[rhs].start < lhs_end && lhs_end < rhs_end)
            {
                overlap[rhs] = max(overlap[rhs], lhs_end);
            } // if
        } // for
    } // for
} // overlap_array


static size_t
cover(size_t const len, size_t const n, struct PMR const pmrs[n], size_t max_cover[len])
{
    size_t inv[len][n + 1];
    inv_array(len, n, pmrs, inv);

    size_t overlap[n];
    overlap_array(n, pmrs, overlap);

    size_t work = 0;

    max_cover[0] = 0;
    for (size_t pos = 1; pos < len; ++pos)
    {
        size_t value = max_cover[pos - 1];
        //work += 1;

        for (size_t i = 0; i < inv[pos][0]; ++i)
        {
            size_t const j = inv[pos][i + 1];
            size_t const count = (pos - pmrs[j].start + 1) / pmrs[j].period;
            size_t const length = pmrs[j].period * count;
            size_t const prev = pos > length ? max_cover[pos - length] : 0;

            value = max(value, prev + length);
            work += 1;

            size_t const overlap_end = min(overlap[j], pos - pmrs[j].period * 2 + 1);
            for (size_t p = pmrs[j].start; p < overlap_end; ++p)
            {
                size_t const count = (pos - p) / pmrs[j].period;
                size_t const length = pmrs[j].period * count;
                size_t const prev = max_cover[pos - length];

                value = max(value, prev + length);
                work += 1;
            } // for
        } // for

        max_cover[pos] = value;
    } // for
    return work;
} // cover


static size_t
cover_q(size_t const len, size_t const n, struct PMR const pmrs[n], size_t max_cover[len])
{
    size_t inv[len][n + 1];
    inv_array(len, n, pmrs, inv);

    size_t work = 0;

    size_t value = 0;
    max_cover[0] = value;
    for (size_t pos = 1; pos < len; ++pos)
    {
        for (size_t i = 0; i < inv[pos][0]; ++i)
        {
            size_t const j = inv[pos][i + 1];

            if (pos - pmrs[j].start + 1 >= 3 * pmrs[j].period)
            {
                size_t const prev = pos > 3 * pmrs[j].period ? pos - 3 * pmrs[j].period : 0;
                value = max(value, max_cover[prev] + 3 * pmrs[j].period);
                work += 1;
            } // if

            size_t const prev = pos > 2 * pmrs[j].period ? pos - 2 * pmrs[j].period : 0;
            value = max(value, max_cover[prev] + 2 * pmrs[j].period);
            work += 1;
        } // for

        max_cover[pos] = value;
    } // for

    return work;
} // cover_q


static size_t
compare_cover(size_t const len, char const word[len])
{
    if (len == 0)
    {
        return 1;
    } // if

    struct PMR pmrs[len];
    size_t const n_pmrs = find_pmrs(len, word, pmrs);
    size_t max_cover[len];
    size_t const work_a = cover(len, n_pmrs, pmrs, max_cover);
    size_t const max_a = max_cover[len - 1];
    size_t const work_b = cover_q(len, n_pmrs, pmrs, max_cover);
    size_t const max_b = max_cover[len - 1];
    if (max_a != max_b)
    {
        printf("%zu %.*s %zu :: %zu vs %zu\n", len, (int) len, word, n_pmrs, max_a, max_b);
    } // if
    return 1;
} // compare_cover


size_t find_runs(char const* word, struct PMR* pmrs);


static size_t g_max_count = 0;
static size_t g_max_depth = 0;
static size_t g_max_sum = 0;


static size_t
count_pmrs(size_t const len, char const word[len])
{
    struct PMR pmrs[len];
    size_t const n = find_pmrs_alt(len, word, pmrs);
    //size_t const n = find_runs(word, pmrs);

    g_max_count = max(g_max_count, n);

    size_t depth[len];
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
    g_max_depth = max(g_max_depth, max_depth);
    g_max_sum = max(g_max_sum, sum);

    return 1;
} // count_pmrs


static size_t
non_isomorphic_n_ary_words(size_t const len,
                           size_t const n,
                           char const alphabet[n],
                           size_t const prev,
                           size_t const current,
                           char word[len],
                           size_t callback(size_t const n, char const word[n]))
{
    if (current >= len)
    {
        return callback(len, word);
    } // if

    if (prev == (size_t) -1)
    {
        word[0] = alphabet[0];
        return non_isomorphic_n_ary_words(len, n, alphabet, 0, 1, word, callback);
    } // if

    size_t count = 0;
    for (size_t i = 0; i <= min(n - 1, prev + 1); ++i)
    {
        word[current] = alphabet[i];
        count += non_isomorphic_n_ary_words(len, n, alphabet, max(i, prev), current + 1, word, callback);
    } // for
    return count;
} // non_isomorphic_n_ary_words


inline static bool
intersect(struct PMR const lhs, struct PMR const rhs)
{
    return rhs.start < lhs.start + lhs.period * lhs.count + lhs.remainder &&
           lhs.start < rhs.start + rhs.period * rhs.count + rhs.remainder;
} // intersect


static size_t
brute_cover(size_t const len,
            char const word[len],
            size_t const n,
            struct PMR const pmrs[n],
            size_t const prev,
            size_t const current,
            size_t const value,
            struct PMR cover[n])
{
    // no more pmrs
    if (current >= n)
    {
        printf("%2zu: [", value);
        for (size_t i = 0; i < prev; ++i)
        {
            printf("(%2zu, %2zu, %2zu), ", cover[i].start, cover[i].period, cover[i].count);
        } // for
        printf("]\n");
        return value;
    } // if

    // skip the current pmr
    size_t local_value = brute_cover(len, word, n, pmrs, prev, current + 1, value, cover);

    // try all instances of the current pmr
    size_t const pmr_end = pmrs[current].start + pmrs[current].period * pmrs[current].count + pmrs[current].remainder;
    for (size_t p = pmrs[current].start; p <= pmr_end - pmrs[current].period * 2; ++p)
    {
        for (size_t c = 2; p + pmrs[current].period * c <= pmr_end; ++c)
        {
            bool safe = true;
            for (size_t i = 0; i < prev; ++i)
            {
                if (intersect(cover[i], (struct PMR) {p, pmrs[current].period, c, 0}))
                {
                    safe = false;
                    break;
                } // if
            } // for

            if (safe)
            {
                cover[prev].start = p;
                cover[prev].period = pmrs[current].period;
                cover[prev].count = c;
                cover[prev].remainder = 0;
                size_t const local = brute_cover(len, word, n, pmrs, prev + 1, current + 1, value + pmrs[current].period * c, cover);
                if (local > local_value)
                {
                    local_value = local;
                } // if
            } // if
        } // for
    } // for

    return local_value;
} // brute_cover


static void
read(FILE* stream)
{
    size_t len = 0;
    char word[WORD_SIZE] = {'\0'};
    size_t n = 0;
    size_t max_cover = 0;

    while (fscanf(stream, "%zu %64s %zu %zu", &len, word, &n, &max_cover) == 4)
    {
        struct PMR pmrs[len];
        size_t const n_pmrs = find_pmrs(len, word, pmrs);
        size_t cover_array[len];
        size_t const work = cover(len, n_pmrs, pmrs, cover_array);
        if (n != n_pmrs || cover_array[len - 1] != max_cover)
        {
            printf("ERROR: %zu %s %zu %zu %zu\n", len, word, n, max_cover, n_pmrs);
            exit(EXIT_FAILURE);
        } // if
    } // while
} // read


static void
read2(FILE* stream)
{
    size_t len = 0;
    char word[24] = {'\0'};
    size_t rrs_len = 0;
    char rrs_word[64] = {'\0'};
    size_t count = 0;

    while (fscanf(stream, "%zu %24s %zu %64s", &len, word, &rrs_len, rrs_word) == 4)
    {
        printf("%zu %s %zu %s\n", len, word, rrs_len, rrs_word);
        count += 1;
    } // while
    fprintf(stderr, "%zu\n", count);
} // read2



int
main(int argc, char* argv[])
{
    read2(stdin);
    return EXIT_SUCCESS;

    /*
    for (size_t i = 30; i <= 30; ++i)
    {
        g_max_count = 0;
        g_max_depth = 0;
        g_max_sum = 0;

        static char word[128] = {'\0'};
        non_isomorphic_n_ary_words(i, 2, "ACGT", -1, 0, word, count_pmrs);
        printf("%zu %zu %zu %zu\n", i, g_max_count, g_max_depth, g_max_sum);
        fflush(stdout);
    } // for
    */
/*
    return -1;

    //printf("%d\n", are_isomorphic_words(4, "ACAB", "XCXC"));
    //return -1;

    //printf("%d\n", intersect((struct PMR) {5, 3, 2, 0}, (struct PMR) {3, 1, 2, 0}));
    //return -1;

    if (argc < 2)
    {
        //fprintf(stderr, "usage: %s length [string]\n", argv[0]);
        read(stdin);
        return EXIT_SUCCESS;
    } // if

    size_t const upto = atoll(argv[1]);

    if (upto == 0)
    {
        if (argc < 3)
        {
            fprintf(stderr, "usage: %s length [string]\n", argv[0]);
            return EXIT_FAILURE;
        } // if

        char const* const word = argv[2];
        size_t const len = strlen(word);
        printf("%zu %s\n", len, word);

        struct PMR pmrs[len];
        size_t const n_pmrs = find_pmrs(len, word, pmrs);
        size_t overlap[len];
        overlap_array(n_pmrs, pmrs, overlap);

        for (size_t i = 0; i < n_pmrs; ++i)
        {
            printf("(%2zu, %2zu, %2zu, %2zu),  # %2zu: %.*s: %2zu\n", pmrs[i].start, pmrs[i].period, pmrs[i].count, pmrs[i].remainder, i, (int) pmrs[i].period, argv[2] + pmrs[i].start, overlap[i]);
        } // for

        for (size_t i = 0; i < len; ++i)
        {
            printf(" %2zu", i);
        } // for
        printf("\n");
        for (size_t i = 0; i < len; ++i)
        {
            printf(" %2c", word[i]);
        } // for
        printf("\n");

        size_t inv[len][n_pmrs + 1];
        size_t const depth = inv_array(len, n_pmrs, pmrs, inv);
        for (size_t i = 1; i <= depth; ++i)
        {
            for (size_t j = 0; j < len; ++j)
            {
                if (inv[j][0] < i)
                {
                    if (i == 1 && inv[j][0] == 0)
                    {
                        printf("  .");
                    } // if
                    else
                    {
                        printf("   ");
                    } // else
                } // if
                else
                {
                    printf(" %2zu", inv[j][i]);
                } // else
            } // for
            printf("\n");
        } // for
        size_t max_cover[len];
        size_t const work = cover(len, n_pmrs, pmrs, max_cover);
        for (size_t i = 0; i < len; ++i)
        {
            printf(" %2zu", max_cover[i]);
        } // for
        printf("\n");
        struct PMR covers[n_pmrs];
        size_t const max_brute = brute_cover(len, word, n_pmrs, pmrs, 0, 0, 0, covers);
        printf("max brute: %zu\n", max_brute);
        printf("work: %zu\n", work);

        return EXIT_SUCCESS;
    } // if

    static char word[WORD_SIZE] = {'\0'};
    for (size_t i = 0; i < (size_t) ((1 << upto) - 1); ++i)
    {
        size_t const len = decode(WORD_SIZE, word, i);
        struct PMR pmrs[len];
        size_t const n_pmrs = find_pmrs(len, word, pmrs);
        size_t cover_array[len];
        size_t const work = cover(len, n_pmrs, pmrs, cover_array);
        printf("%zu %s %zu %zu %zu\n", len, word, n_pmrs, cover_array[len - 1], work);
    } // for
*/
    return EXIT_SUCCESS;
} // main
