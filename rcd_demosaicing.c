/**
* RATIO CORRECTED DEMOSAICING
* Luis Sanz Rodríguez (luis.sanz.rodriguez(at)gmail(dot)com)
*
* Release 2.3 @ 171125
*/

void CLASS rcd_demosaicing()
{
    int row, col, indx, c;
    int w1 = width, w2 = 2 * width, w3 = 3 * width, w4 = 4 * width;
    float (*cfa), (*rgb)[3];

    //Tolerance to avoid dividing by zero
    static const float eps = 1e-5, epssq = 1e-10;

    //Gradients
    float N_Grad, E_Grad, W_Grad, S_Grad, NW_Grad, NE_Grad, SW_Grad, SE_Grad;

    //Pixel estimation
    float N_Est, E_Est, W_Est, S_Est, NW_Est, NE_Est, SW_Est, SE_Est, V_Est, H_Est, P_Est, Q_Est;

    //Directional discrimination
    float V_Stat, H_Stat, P_Stat, Q_Stat;
    float VH_Central_Value, VH_Neighbourhood_Value, PQ_Central_Value, PQ_Neighbourhood_Value;
    float ( *VH_Dir ), VH_Disc, ( *PQ_Dir ), PQ_Disc;

    //Low pass filter
    float ( *lpf );


    //Convert the CFA to floating point buffers
    cfa = ( float ( * )    ) calloc( width * height, sizeof *cfa ); merror ( cfa, "rcd_demosaicing_171125()" );
    rgb = ( float ( * )[3] ) calloc( width * height, sizeof *rgb ); merror ( rgb, "rcd_demosaicing_171125()" );

    if ( verbose ) fprintf ( stderr, _( "RCD interpolation...\n" ) );

    for ( row = 0; row < height; row++ ) {
        for ( col = 0, indx = row * width + col; col < width; col++, indx++ ) {

            cfa[indx] = rgb[indx][FC( row, col )] = (float)image[indx][FC( row, col )] / 65535.f;

        }
    }


    /**
    * STEP 1: Find vertical and horizontal interpolation directions
    */

    VH_Dir = ( float ( * ) ) calloc( width * height, sizeof *VH_Dir ); merror ( VH_Dir, "rcd_demosaicing_171125()" );

    // Step 1.1: Calculate vertical and horizontal local discrimination
    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4, indx = row * width + col; col < width - 4; col++, indx++ ) {

            V_Stat = MAX( - 18.f  *  cfa[indx] * cfa[indx - w1] - 18.f * cfa[indx] * cfa[indx + w1] - 36.f * cfa[indx] * cfa[indx - w2] - 36.f * cfa[indx] * cfa[indx + w2] + 18.f * cfa[indx] * cfa[indx - w3] + 18.f * cfa[indx] * cfa[indx + w3] - 2.f * cfa[indx] * cfa[indx - w4] - 2.f * cfa[indx] * cfa[indx + w4] + 38.f * cfa[indx] * cfa[indx] - 70.f * cfa[indx - w1] * cfa[indx + w1] - 12.f * cfa[indx - w1] * cfa[indx - w2] + 24.f * cfa[indx - w1] * cfa[indx + w2] - 38.f * cfa[indx - w1] * cfa[indx - w3] + 16.f * cfa[indx - w1] * cfa[indx + w3] + 12.f * cfa[indx - w1] * cfa[indx - w4] - 6.f * cfa[indx - w1] * cfa[indx + w4] + 46.f * cfa[indx - w1] * cfa[indx - w1] + 24.f * cfa[indx + w1] * cfa[indx - w2] - 12.f * cfa[indx + w1] * cfa[indx + w2] + 16.f * cfa[indx + w1] * cfa[indx - w3] - 38.f * cfa[indx + w1] * cfa[indx + w3] - 6.f * cfa[indx + w1] * cfa[indx - w4] + 12.f * cfa[indx + w1] * cfa[indx + w4] + 46.f * cfa[indx + w1] * cfa[indx + w1] + 14.f * cfa[indx - w2] * cfa[indx + w2] - 12.f * cfa[indx - w2] * cfa[indx + w3] - 2.f * cfa[indx - w2] * cfa[indx - w4] + 2.f * cfa[indx - w2] * cfa[indx + w4] + 11.f * cfa[indx - w2] * cfa[indx - w2] - 12.f * cfa[indx + w2] * cfa[indx - w3] + 2.f * cfa[indx + w2] * cfa[indx - w4] - 2.f * cfa[indx + w2] * cfa[indx + w4] + 11.f * cfa[indx + w2] * cfa[indx + w2] + 2.f * cfa[indx - w3] * cfa[indx + w3] - 6.f * cfa[indx - w3] * cfa[indx - w4] + 10.f * cfa[indx - w3] * cfa[indx - w3] - 6.f * cfa[indx + w3] * cfa[indx + w4] + 10.f * cfa[indx + w3] * cfa[indx + w3] + 1.f * cfa[indx - w4] * cfa[indx - w4] + 1.f * cfa[indx + w4] * cfa[indx + w4], epssq );
            H_Stat = MAX( - 18.f  *  cfa[indx] * cfa[indx -  1] - 18.f * cfa[indx] * cfa[indx +  1] - 36.f * cfa[indx] * cfa[indx -  2] - 36.f * cfa[indx] * cfa[indx +  2] + 18.f * cfa[indx] * cfa[indx -  3] + 18.f * cfa[indx] * cfa[indx +  3] - 2.f * cfa[indx] * cfa[indx -  4] - 2.f * cfa[indx] * cfa[indx +  4] + 38.f * cfa[indx] * cfa[indx] - 70.f * cfa[indx -  1] * cfa[indx +  1] - 12.f * cfa[indx -  1] * cfa[indx -  2] + 24.f * cfa[indx -  1] * cfa[indx +  2] - 38.f * cfa[indx -  1] * cfa[indx -  3] + 16.f * cfa[indx -  1] * cfa[indx +  3] + 12.f * cfa[indx -  1] * cfa[indx -  4] - 6.f * cfa[indx -  1] * cfa[indx +  4] + 46.f * cfa[indx -  1] * cfa[indx -  1] + 24.f * cfa[indx +  1] * cfa[indx -  2] - 12.f * cfa[indx +  1] * cfa[indx +  2] + 16.f * cfa[indx +  1] * cfa[indx -  3] - 38.f * cfa[indx +  1] * cfa[indx +  3] - 6.f * cfa[indx +  1] * cfa[indx -  4] + 12.f * cfa[indx +  1] * cfa[indx +  4] + 46.f * cfa[indx +  1] * cfa[indx +  1] + 14.f * cfa[indx -  2] * cfa[indx +  2] - 12.f * cfa[indx -  2] * cfa[indx +  3] - 2.f * cfa[indx -  2] * cfa[indx -  4] + 2.f * cfa[indx -  2] * cfa[indx +  4] + 11.f * cfa[indx -  2] * cfa[indx -  2] - 12.f * cfa[indx +  2] * cfa[indx -  3] + 2.f * cfa[indx +  2] * cfa[indx -  4] - 2.f * cfa[indx +  2] * cfa[indx +  4] + 11.f * cfa[indx +  2] * cfa[indx +  2] + 2.f * cfa[indx -  3] * cfa[indx +  3] - 6.f * cfa[indx -  3] * cfa[indx -  4] + 10.f * cfa[indx -  3] * cfa[indx -  3] - 6.f * cfa[indx +  3] * cfa[indx +  4] + 10.f * cfa[indx +  3] * cfa[indx +  3] + 1.f * cfa[indx -  4] * cfa[indx -  4] + 1.f * cfa[indx +  4] * cfa[indx +  4], epssq );

            VH_Dir[indx] = V_Stat / ( V_Stat + H_Stat );

        }
    }


    /**
    * STEP 2: Calculate the low pass filter
    */

    // Step 2.1: Low pass filter incorporating green, red and blue local samples from the raw data
    lpf = ( float ( * ) ) calloc( width * height, sizeof *lpf ); merror ( lpf, "rcd_demosaicing_171125()" );

    for ( row = 2; row < height - 2; row++ ) {
        for ( col = 2 + FC( row, 0 ) & 1, indx = row * width + col; col < width - 2; col += 2, indx += 2 ) {

            lpf[indx] = 0.25f * cfa[indx] + 0.125f * ( cfa[indx - w1] + cfa[indx + w1] + cfa[indx - 1] + cfa[indx + 1] ) + 0.0625f * ( cfa[indx - w1 - 1] + cfa[indx - w1 + 1] + cfa[indx + w1 - 1] + cfa[indx + w1 + 1] );

        }
    }


    /**
    * STEP 3: Populate the green channel
    */

    // Step 3.1: Populate the green channel at blue and red CFA positions
    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + FC( row, 0 ) & 1, indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            // Refined vertical and horizontal local discrimination
            VH_Central_Value   = VH_Dir[indx];
            VH_Neighbourhood_Value = 0.25f * ( VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1] + VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1] );

            VH_Disc = ( fabs( 0.5f - VH_Central_Value ) < fabs( 0.5f - VH_Neighbourhood_Value ) ) ? VH_Neighbourhood_Value : VH_Central_Value;

            // Cardinal gradients
            N_Grad = eps + fabs( cfa[indx - w1] - cfa[indx + w1] ) + fabs( cfa[indx] - cfa[indx - w2] ) + fabs( cfa[indx - w1] - cfa[indx - w3] ) + fabs( cfa[indx - w2] - cfa[indx - w4] );
            S_Grad = eps + fabs( cfa[indx + w1] - cfa[indx - w1] ) + fabs( cfa[indx] - cfa[indx + w2] ) + fabs( cfa[indx + w1] - cfa[indx + w3] ) + fabs( cfa[indx + w2] - cfa[indx + w4] );
            W_Grad = eps + fabs( cfa[indx -  1] - cfa[indx +  1] ) + fabs( cfa[indx] - cfa[indx -  2] ) + fabs( cfa[indx -  1] - cfa[indx -  3] ) + fabs( cfa[indx -  2] - cfa[indx -  4] );
            E_Grad = eps + fabs( cfa[indx +  1] - cfa[indx -  1] ) + fabs( cfa[indx] - cfa[indx +  2] ) + fabs( cfa[indx +  1] - cfa[indx +  3] ) + fabs( cfa[indx +  2] - cfa[indx +  4] );

            // Cardinal pixel estimations
            N_Est = cfa[indx - w1] * ( 1.f + ( lpf[indx] - lpf[indx - w2] ) / ( eps + lpf[indx] + lpf[indx - w2] ) );
            S_Est = cfa[indx + w1] * ( 1.f + ( lpf[indx] - lpf[indx + w2] ) / ( eps + lpf[indx] + lpf[indx + w2] ) );
            W_Est = cfa[indx -  1] * ( 1.f + ( lpf[indx] - lpf[indx -  2] ) / ( eps + lpf[indx] + lpf[indx -  2] ) );
            E_Est = cfa[indx +  1] * ( 1.f + ( lpf[indx] - lpf[indx +  2] ) / ( eps + lpf[indx] + lpf[indx +  2] ) );

            // Vertical and horizontal estimations
            V_Est = ( S_Grad * N_Est + N_Grad * S_Est ) / ( N_Grad + S_Grad );
            H_Est = ( W_Grad * E_Est + E_Grad * W_Est ) / ( E_Grad + W_Grad );

            // G@B and G@R interpolation
            rgb[indx][1] = LIM( VH_Disc * H_Est + ( 1.f - VH_Disc ) * V_Est, 0.f, 1.f );

        }
    }

    free( lpf );


    /**
    * STEP 4: Populate the red and blue channels
    */

    // Step 4.1: Calculate P/Q diagonal local discrimination
    PQ_Dir = ( float ( * ) ) calloc( width * height, sizeof *PQ_Dir ); merror ( PQ_Dir, "rcd_demosaicing_171125()" );

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + FC( row, 0 ) & 1, indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            P_Stat = MAX( - 18.f * cfa[indx] * cfa[indx - w1 - 1] - 18.f * cfa[indx] * cfa[indx + w1 + 1] - 36.f * cfa[indx] * cfa[indx - w2 - 2] - 36.f * cfa[indx] * cfa[indx + w2 + 2] + 18.f * cfa[indx] * cfa[indx - w3 - 3] + 18.f * cfa[indx] * cfa[indx + w3 + 3] - 2.f * cfa[indx] * cfa[indx - w4 - 4] - 2.f * cfa[indx] * cfa[indx + w4 + 4] + 38.f * cfa[indx] * cfa[indx] - 70.f * cfa[indx - w1 - 1] * cfa[indx + w1 + 1] - 12.f * cfa[indx - w1 - 1] * cfa[indx - w2 - 2] + 24.f * cfa[indx - w1 - 1] * cfa[indx + w2 + 2] - 38.f * cfa[indx - w1 - 1] * cfa[indx - w3 - 3] + 16.f * cfa[indx - w1 - 1] * cfa[indx + w3 + 3] + 12.f * cfa[indx - w1 - 1] * cfa[indx - w4 - 4] - 6.f * cfa[indx - w1 - 1] * cfa[indx + w4 + 4] + 46.f * cfa[indx - w1 - 1] * cfa[indx - w1 - 1] + 24.f * cfa[indx + w1 + 1] * cfa[indx - w2 - 2] - 12.f * cfa[indx + w1 + 1] * cfa[indx + w2 + 2] + 16.f * cfa[indx + w1 + 1] * cfa[indx - w3 - 3] - 38.f * cfa[indx + w1 + 1] * cfa[indx + w3 + 3] - 6.f * cfa[indx + w1 + 1] * cfa[indx - w4 - 4] + 12.f * cfa[indx + w1 + 1] * cfa[indx + w4 + 4] + 46.f * cfa[indx + w1 + 1] * cfa[indx + w1 + 1] + 14.f * cfa[indx - w2 - 2] * cfa[indx + w2 + 2] - 12.f * cfa[indx - w2 - 2] * cfa[indx + w3 + 3] - 2.f * cfa[indx - w2 - 2] * cfa[indx - w4 - 4] + 2.f * cfa[indx - w2 - 2] * cfa[indx + w4 + 4] + 11.f * cfa[indx - w2 - 2] * cfa[indx - w2 - 2] - 12.f * cfa[indx + w2 + 2] * cfa[indx - w3 - 3] + 2 * cfa[indx + w2 + 2] * cfa[indx - w4 - 4] - 2.f * cfa[indx + w2 + 2] * cfa[indx + w4 + 4] + 11.f * cfa[indx + w2 + 2] * cfa[indx + w2 + 2] + 2.f * cfa[indx - w3 - 3] * cfa[indx + w3 + 3] - 6.f * cfa[indx - w3 - 3] * cfa[indx - w4 - 4] + 10.f * cfa[indx - w3 - 3] * cfa[indx - w3 - 3] - 6.f * cfa[indx + w3 + 3] * cfa[indx + w4 + 4] + 10.f * cfa[indx + w3 + 3] * cfa[indx + w3 + 3] + 1.f * cfa[indx - w4 - 4] * cfa[indx - w4 - 4] + 1.f * cfa[indx + w4 + 4] * cfa[indx + w4 + 4], epssq );
            Q_Stat = MAX( - 18.f * cfa[indx] * cfa[indx + w1 - 1] - 18.f * cfa[indx] * cfa[indx - w1 + 1] - 36.f * cfa[indx] * cfa[indx + w2 - 2] - 36.f * cfa[indx] * cfa[indx - w2 + 2] + 18.f * cfa[indx] * cfa[indx + w3 - 3] + 18.f * cfa[indx] * cfa[indx - w3 + 3] - 2.f * cfa[indx] * cfa[indx + w4 - 4] - 2.f * cfa[indx] * cfa[indx - w4 + 4] + 38.f * cfa[indx] * cfa[indx] - 70.f * cfa[indx + w1 - 1] * cfa[indx - w1 + 1] - 12.f * cfa[indx + w1 - 1] * cfa[indx + w2 - 2] + 24.f * cfa[indx + w1 - 1] * cfa[indx - w2 + 2] - 38.f * cfa[indx + w1 - 1] * cfa[indx + w3 - 3] + 16.f * cfa[indx + w1 - 1] * cfa[indx - w3 + 3] + 12.f * cfa[indx + w1 - 1] * cfa[indx + w4 - 4] - 6.f * cfa[indx + w1 - 1] * cfa[indx - w4 + 4] + 46.f * cfa[indx + w1 - 1] * cfa[indx + w1 - 1] + 24.f * cfa[indx - w1 + 1] * cfa[indx + w2 - 2] - 12.f * cfa[indx - w1 + 1] * cfa[indx - w2 + 2] + 16.f * cfa[indx - w1 + 1] * cfa[indx + w3 - 3] - 38.f * cfa[indx - w1 + 1] * cfa[indx - w3 + 3] - 6.f * cfa[indx - w1 + 1] * cfa[indx + w4 - 4] + 12.f * cfa[indx - w1 + 1] * cfa[indx - w4 + 4] + 46.f * cfa[indx - w1 + 1] * cfa[indx - w1 + 1] + 14.f * cfa[indx + w2 - 2] * cfa[indx - w2 + 2] - 12.f * cfa[indx + w2 - 2] * cfa[indx - w3 + 3] - 2.f * cfa[indx + w2 - 2] * cfa[indx + w4 - 4] + 2.f * cfa[indx + w2 - 2] * cfa[indx - w4 + 4] + 11.f * cfa[indx + w2 - 2] * cfa[indx + w2 - 2] - 12.f * cfa[indx - w2 + 2] * cfa[indx + w3 - 3] + 2 * cfa[indx - w2 + 2] * cfa[indx + w4 - 4] - 2.f * cfa[indx - w2 + 2] * cfa[indx - w4 + 4] + 11.f * cfa[indx - w2 + 2] * cfa[indx - w2 + 2] + 2.f * cfa[indx + w3 - 3] * cfa[indx - w3 + 3] - 6.f * cfa[indx + w3 - 3] * cfa[indx + w4 - 4] + 10.f * cfa[indx + w3 - 3] * cfa[indx + w3 - 3] - 6.f * cfa[indx - w3 + 3] * cfa[indx - w4 + 4] + 10.f * cfa[indx - w3 + 3] * cfa[indx - w3 + 3] + 1.f * cfa[indx + w4 - 4] * cfa[indx + w4 - 4] + 1.f * cfa[indx - w4 + 4] * cfa[indx - w4 + 4], epssq );

            PQ_Dir[indx] = P_Stat / ( P_Stat + Q_Stat );

        }
    }

    // Step 4.2: Populate the red and blue channels at blue and red CFA positions
    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + FC( row, 0 ) & 1, indx = row * width + col, c = 2 - FC( row, col ); col < width - 4; col += 2, indx += 2 ) {

            // Refined P/Q diagonal local discrimination
            PQ_Central_Value   = PQ_Dir[indx];
            PQ_Neighbourhood_Value = 0.25f * ( PQ_Dir[indx - w1 - 1] + PQ_Dir[indx - w1 + 1] + PQ_Dir[indx + w1 - 1] + PQ_Dir[indx + w1 + 1] );

            PQ_Disc = ( fabs( 0.5f - PQ_Central_Value ) < fabs( 0.5f - PQ_Neighbourhood_Value ) ) ? PQ_Neighbourhood_Value : PQ_Central_Value;

            // Diagonal gradients
            NW_Grad = eps + fabs( rgb[indx - w1 - 1][c] - rgb[indx + w1 + 1][c] ) + fabs( rgb[indx - w1 - 1][c] - rgb[indx - w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 - 2][1] );
            NE_Grad = eps + fabs( rgb[indx - w1 + 1][c] - rgb[indx + w1 - 1][c] ) + fabs( rgb[indx - w1 + 1][c] - rgb[indx - w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 + 2][1] );
            SW_Grad = eps + fabs( rgb[indx + w1 - 1][c] - rgb[indx - w1 + 1][c] ) + fabs( rgb[indx + w1 - 1][c] - rgb[indx + w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 - 2][1] );
            SE_Grad = eps + fabs( rgb[indx + w1 + 1][c] - rgb[indx - w1 - 1][c] ) + fabs( rgb[indx + w1 + 1][c] - rgb[indx + w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 + 2][1] );

            // Diagonal colour differences
            NW_Est = rgb[indx - w1 - 1][c] - rgb[indx - w1 - 1][1];
            NE_Est = rgb[indx - w1 + 1][c] - rgb[indx - w1 + 1][1];
            SW_Est = rgb[indx + w1 - 1][c] - rgb[indx + w1 - 1][1];
            SE_Est = rgb[indx + w1 + 1][c] - rgb[indx + w1 + 1][1];

            // P/Q estimations
            P_Est = ( NW_Grad * SE_Est + SE_Grad * NW_Est ) / ( NW_Grad + SE_Grad );
            Q_Est = ( NE_Grad * SW_Est + SW_Grad * NE_Est ) / ( NE_Grad + SW_Grad );

            // R@B and B@R interpolation
            rgb[indx][c] = LIM( rgb[indx][1] + ( 1.f - PQ_Disc ) * P_Est + PQ_Disc * Q_Est, 0.f, 1.f );

        }
    }

    free( PQ_Dir );

    // Step 4.3: Populate the red and blue channels at green CFA positions
    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + FC( row, 1 ) & 1, indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            // Refined vertical and horizontal local discrimination
            VH_Central_Value   = VH_Dir[indx];
            VH_Neighbourhood_Value = 0.25f * ( VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1] + VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1] );

            VH_Disc = ( fabs( 0.5f - VH_Central_Value ) < fabs( 0.5f - VH_Neighbourhood_Value ) ) ? VH_Neighbourhood_Value : VH_Central_Value;

            for ( c = 0; c <= 2; c += 2 ) {

                // Cardinal gradients
                N_Grad = eps + fabs( rgb[indx][1] - rgb[indx - w2][1] ) + fabs( rgb[indx - w1][c] - rgb[indx + w1][c] ) + fabs( rgb[indx - w1][c] - rgb[indx - w3][c] );
                S_Grad = eps + fabs( rgb[indx][1] - rgb[indx + w2][1] ) + fabs( rgb[indx + w1][c] - rgb[indx - w1][c] ) + fabs( rgb[indx + w1][c] - rgb[indx + w3][c] );
                W_Grad = eps + fabs( rgb[indx][1] - rgb[indx -  2][1] ) + fabs( rgb[indx -  1][c] - rgb[indx +  1][c] ) + fabs( rgb[indx -  1][c] - rgb[indx -  3][c] );
                E_Grad = eps + fabs( rgb[indx][1] - rgb[indx +  2][1] ) + fabs( rgb[indx +  1][c] - rgb[indx -  1][c] ) + fabs( rgb[indx +  1][c] - rgb[indx +  3][c] );

                // Cardinal colour differences
                N_Est = rgb[indx - w1][c] - rgb[indx - w1][1];
                S_Est = rgb[indx + w1][c] - rgb[indx + w1][1];
                W_Est = rgb[indx -  1][c] - rgb[indx -  1][1];
                E_Est = rgb[indx +  1][c] - rgb[indx +  1][1];

                // Vertical and horizontal estimations
                V_Est = ( N_Grad * S_Est + S_Grad * N_Est ) / ( N_Grad + S_Grad );
                H_Est = ( E_Grad * W_Est + W_Grad * E_Est ) / ( E_Grad + W_Grad );

                // R@G and B@G interpolation
                rgb[indx][c] = LIM( rgb[indx][1] + ( 1.f - VH_Disc ) * V_Est + VH_Disc * H_Est, 0.f, 1.f );

            }
        }
    }


    //Revert back the floating point buffers to dcraw (int)image
    for ( row = 0; row < height; row++ ) {
        for ( col = 0, indx = row * width + col; col < width; col++, indx++ ) {

            FORC3 image[indx][c] = (ushort)CLIP( 65535.f * rgb[indx][c] );

        }
    }

    border_interpolate( 4 );

    free( VH_Dir );

}
