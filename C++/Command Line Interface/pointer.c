#include <stdio.h>

// Delete after fixing call_by_reference
int x;

/**
 * @brief Demonstrates basic usage of pointers to access a local variable.
 *
 * This function takes an integer by value, obtains its address,
 * and prints both the value and the address in two ways: directly
 * and via the pointer.
 *
 * @param x Value to demonstrate pointer basics on.
 */
void basic_pointer (int x)
{
    // Pointer to hold the address of x
    int* address_of_x;

    // Store the address of x in the pointer
    address_of_x = &x;

    // Print the value of x directly
    printf("The value of x is %d\n", x /* TODO */);
    // Print the address of x (directly using &x)
    printf("The address of x is %p\n", (void*)&x /* TODO */);
    // Print the address of x (via the pointer)
    printf("The address of x via address_of_x is %p\n", (void*)address_of_x /* TODO */);
    // Print the value of x by dereferencing the pointer
    printf("The value of x via address_of_x is %d\n", *address_of_x /* TODO */);
}

/**
 * @brief Shows how assignment via pointer works with two variables.
 *
 * Copies the value of x into y, then modifies x and shows how y can
 * be re-assigned via the pointer to x.
 *
 * @param x Initial value to work with.
 */
void basic_pointer2 (int x)
{
    // Initialize pointer to the address of x
    int* address_of_x = &x /* TODO */;
    // Another variable y gets the value of x
    int y = x /* TODO */;

    // Print the current value of y
    printf("The value of y is %d\n", y /* TODO */);

    // Change the value of x directly
    x = 10;
    // Assign y from the new value of x via the pointer
    y = *address_of_x;

    // Print the updated value of y
    printf("The value of y is %d\n", y /* TODO */);
}

/**
 * @brief Modifies the value of its local parameter via pointer.
 *
 * Demonstrates writing through a pointer to change the underlying
 * integer value, then prints after each modification.
 *
 * @param x Value to be modified via its address.
 */
void basic_pointer_changeValue (int x)
{
    // Pointer to x’s address
    int* address_of_x = &x /* TODO */;

    // Change the value of x to 10 via the pointer
    *address_of_x = 10;
    printf("x = %d\n", x /* TODO */);

    // Change the value of x again to 20 via the pointer
    *address_of_x = 20;
    printf("x = %d\n", x);
}

/**
 * @brief Changes the caller's variable by reference.
 *
 * This function receives a pointer to an integer and sets the
 * pointed-to value to 200, demonstrating “call by reference.”
 *
 * @param x Pointer to an integer to be modified.
 */
void call_by_reference (int* x)
{
    // Change the value stored at the address in x to 200
    *x = 200;
}

/**
 * @brief Entry point of the program.
 *
 * Calls each demonstration function in turn to show pointer usage
 * and reference-based modification.
 *
 * @return Exit status code (0 = success).
 */
int main (void)
{
    // Initialize x in main
    int x = 5;

    // Demonstrate basic pointer operations
    basic_pointer(x);
    // Demonstrate assignment via pointer
    basic_pointer2(x);
    // Demonstrate modifying a local copy via pointer
    basic_pointer_changeValue(x);

    // Show x before call-by-reference
    printf("The value of x before call_by_reference is %d\n", x);
    // Modify x in main via call_by_reference
    call_by_reference(&x);
    // Show x after call-by-reference
    printf("The value of x after call_by_reference is %d\n", x);

    return 0;
}
