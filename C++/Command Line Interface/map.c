/*=============================================================================
 * File:        map.c
 * Author:      Lawrence‑Alagappan‑Somasundaram
 * Purpose:     Demonstrate how to store a **single, outward‑facing** compass
 *              direction in each cell of a 3 × 3 grid by means of bit‑flags and
 *              to print the grid in a readable form.
 *
 * Grid layout (row = x, column = y):
 *
 *      (0,0) NW   (0,1) N   (0,2) NE
 *      (1,0) W    (1,1) ·   (1,2) E
 *      (2,0) SW   (2,1) S   (2,2) SE
 *
 *   • Every cell except the centre (1,1) has exactly **one** valid direction:
 *       – Edge‑centre cells store a cardinal direction (N, E, S, W).
 *       – Corner cells store the matching diagonal (NW, NE, SW, SE).
 *   • Requests that are out‑of‑bounds, contain no bits, contain the wrong bits
 *     for the addressed cell or contain too many/opposite bits are **ignored**.
 *
 * Building & running (POSIX shell):
 *     $ gcc map.c -o map && ./map
 *
 *   The `main()` function prints the grid twice: the first time with uniform
 *   spacing, the second time with a staggered layout to showcase formatting.
 *============================================================================*/

 #include <stdio.h>    /* putchar, fputs */

/*==========================================================================*/
/*  Cardinal‑direction bit‑flags                                            */
/*==========================================================================*/

/**
 * @enum cardd
 * @brief Single‑bit masks for the four primary compass points.
 *
 * Bits (LSB..MSB): N E S W → 0 1 2 3
 *
 * Diagonals are formed at **run‑time** by OR‑ing two appropriate bits; they are
 * *not* members of the enum so that the enumeration remains small and every
 * diagonal is guaranteed to be a unique combination.
 */
 
 typedef enum cardd {
     N = 1 << 0,  /* 0001 */
     E = 1 << 1,  /* 0010 */
     S = 1 << 2,  /* 0100 */
     W = 1 << 3   /* 1000 */
 } cardd;
 
/*--------------------------------------------------------------------------*/
/*  Global 3 × 3 direction grid                                             */
/*--------------------------------------------------------------------------*/

/**
 * @var map
 * @brief Global 3 × 3 array that stores the direction mask for each cell.
 *
 * Initialised to 0 (no direction).  Index order is **row‑major** `[x][y]` to
 * match the *(row,column)* layout used throughout the code.
 */
 static cardd map[3][3] = { { 0 } };
 
 /*==========================================================================*/
/*  set_dir function                                                        */
/*==========================================================================*/

/**
 * @brief Validate and store a direction in cell *(x, y)*.
 *
 * @param x   Row index `0‥2` (top = 0, bottom = 2).
 * @param y   Column index `0‥2` (left = 0, right = 2).
 * @param dir Bit‑flag direction (combination of `N`, `E`, `S`, `W`).
 *
 * The function computes the **single** correct bit‑pattern for coordinates
 * *(x, y)* and writes `dir` into `map[x][y]` **only** when the arguments are
 * valid and `dir` matches that pattern.  Otherwise the function returns
 * without side effects.
 *
 * @note The centre cell (1,1) never accepts a value because its expected mask
 *       is 0.
 */
 void set_dir(int x, int y, enum cardd dir)
 {
     /* 1. quick rejects ---------------------------------------------------- */
     if (x < 0 || x > 2 || y < 0 || y > 2) return; /* out of range */
     if (dir == 0)                               return; /* nothing to set */
 
     /* 2. derive the single valid bit‑pattern for this cell ---------------- */
     enum cardd expect = 0;
 
     /* vertical component */
     if      (x < 1) expect |= N; /* top row    */
     else if (x > 1) expect |= S; /* bottom row */
     /* middle row contributes no vertical bit */
 
     /* horizontal component */
     if      (y < 1) expect |= W; /* left column */
     else if (y > 1) expect |= E; /* right column */
 
     /* The centre cell (1,1) yields expect == 0 – it never stores anything. */
     if (expect == 0) return;
 
     /* 3. accept only the correct value ------------------------------------ */
     if (dir == expect) {
         map[x][y] = dir;
     }
 }
 
 /*==========================================================================*/
/*  set_dir function                                                        */
/*==========================================================================*/
  /**
 * @brief Print the current 3 × 3 map in a human‑readable form.
 *
 * The function translates each direction mask into its two‑letter code and
 * prints a blank cell as `'0'`.  A static call counter alternates the number
 * of spaces between columns so the two calls in `main()` visually differ:
 *
 *   • Call #1 → three spaces between **all** cells.
 *   • Call ≥2 → even rows get two spaces, odd rows keep three.
 *
 * @note The switch expression is explicitly cast to `int` to avoid the
 *       `‑Wswitch‑enum` warning that would otherwise arise from using diagonal
 *       masks (which are not part of `enum cardd`).  Functionality is unchanged.
 */
 void show_map(void)
 {
     static int call = 0;
     ++call;
 
     for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            switch ((int)map[i][j]) {     
                case N:     fputs("N",  stdout); break;
                case E:     fputs("E",  stdout); break;
                case S:     fputs("S",  stdout); break;
                case W:     fputs("W",  stdout); break;
                case N|E:   fputs("NE", stdout); break;
                case N|W:   fputs("NW", stdout); break;
                case S|E:   fputs("SE", stdout); break;
                case S|W:   fputs("SW", stdout); break;
                default:    fputs("0",  stdout); break;
            }
 
             /* Construct the horizontal spacing between columns */
             if (j < 2) {
                 int pad = (call == 1) ? 3           /* first map */
                                       : ((i % 2) ? 3 : 2); /* thereafter */
                 while (pad--) putchar(' ');
             }
         }
         putchar('\n');
     }
     putchar('\n');
 }
 
 /*--------------------------------------------------------------------------
  * Function: main – DO NOT MODIFY.
  *--------------------------------------------------------------------------*/
 int main (void)
 {
     // You are not allowed to change anything in this function!
     set_dir(0, 1, N);
     set_dir(1, 0, W);
     set_dir(1, 4, W);
     set_dir(1, 2, E);
     set_dir(2, 1, S);
 
     show_map();
 
     set_dir(0, 0, N|W);
     set_dir(0, 2, N|E);
     set_dir(0, 2, N|S);
     set_dir(2, 0, S|W);
     set_dir(2, 2, S|E);
     set_dir(2, 2, E|W);
     set_dir(1, 3, N|S|E);
     set_dir(1, 1, N|S|E|W);
 
     show_map();
 
     return 0;
 }
 