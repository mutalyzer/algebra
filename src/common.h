#ifndef GVA_COMMON_H
#define GVA_COMMON_H


// Generic abs function. Evaluates its parameters multiple times.
#define ABS(x) (((x) >= 0) ? (x) : -(x))


// Generic max function. Evaluates its parameters multiple times.
#define MAX(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))


// Generic min function. Evaluates its parameters multiple times.
#define MIN(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))


#endif // GVA_COMMON_H
