#include <stdio.h>
#include <stdlib.h>

void c_convolve(const int nsources,             /*scalar - number of sources*/
                const int noutlets,             /*scalar - length of subset*/
                const int subset_length,        /*scalar - length of subset*/
                const int x_size,
                const int* source2outlet_ind,   /*1d array - source to outlet mapping*/
                const int* source_y_ind,        /*1d array - source y location*/
                const int* source_x_ind,        /*1d array - source x location*/
                const int* source_time_offset,  /*1d array - source time offset*/
                const double* unit_hydrograph,  /*2d array[times][sources] - unit hydrographs*/
                const double* aggrunin,         /*2d array[ysize][xsize] - vic runoff flux*/
                double* ring)                   /*2d array[times][outlets] - convolution ring*/
{
    const double* runin;                      /*pointer to sources runoff flux*/
    int s, i, j;                              /*counters*/
    int y, x, offset, outlet;                 /*2d indicies*/
    int xyind, rind, uhind;                   /*1d indicies*/

    /*Loop through all sources*/
    for (s = 0; s < nsources; s++) {

        outlet = source2outlet_ind[s];
        y = source_y_ind[s];
        x = source_x_ind[s];
        offset = source_time_offset[s];

        //1d index location
        //2d-->1d indexing goes like this:  ind = y*x_size + x
        xyind = y*x_size + x;

        runin = &aggrunin[xyind];

        // if (*runin > 10000.0) {
        //     printf("runin is big: %f\n", *runin);
        //     printf("xyind: %i\n", xyind);
        //     printf("y: %i\n", y);
        //     printf("x: %i\n", x);

        /* Do the convolution */
        for (i = 0; i < subset_length; i++) {
            j = i + offset;

            //1d index locations
            rind = j * noutlets + outlet;
            uhind = i * nsources + s;

            ring[rind] += unit_hydrograph[uhind] * *runin;
        }
    }
}
