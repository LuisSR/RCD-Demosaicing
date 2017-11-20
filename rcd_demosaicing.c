/**
* RATIO CORRECTED DEMOSAICING
* 2011-2017 Luis Sanz Rodríguez (luis.sanz.rodriguez(at)gmail(dot)com)
*
* Release 2.2 @ 171117
*/

void CLASS rcd_demosaicing_171117()
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
    float VH_Central_Value, VH_Neighbour_Value, PQ_Central_Value, PQ_Neighbour_Value;
    float ( *VH_Dir ), ( *VH_Disc ), ( *PQ_Dir ), ( *PQ_Disc );

    //Low pass filter
    float ( *lpf );

    //Convert the CFA to floating point buffers
    cfa = ( float ( * ) ) calloc( width * height, sizeof *cfa ); merror ( cfa, "rcd_demosaicing_171117()" );
    rgb = ( float ( * )[3] ) calloc( width * height, sizeof *rgb ); merror ( rgb, "rcd_demosaicing_171117()" );

    if ( verbose ) fprintf ( stderr, _( "RCD interpolation...\n" ) );

    for ( row = 0; row < height; row++ ) {
        for ( col = 0, indx = row * width + col; col < width; col++, indx++ ) {

            cfa[indx] = rgb[indx][FC(row,col)] = (float) image[indx][FC(row,col)] / 65535.f;

        }
    }


    /**
    * STEP 1: Find cardinal and diagonal interpolation directions
    */

    VH_Dir = ( float ( * ) ) calloc( width * height, sizeof *VH_Dir ); merror ( VH_Dir, "rcd_demosaicing_171117()" );
    PQ_Dir = ( float ( * ) ) calloc( width * height, sizeof *PQ_Dir ); merror ( PQ_Dir, "rcd_demosaicing_171117()" );

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4, indx = row * width + col; col < width - 4; col++, indx++ ) {

            //Calculate V/H local discrimination
            V_Stat = epssq - 18.0f  *  cfa[indx] * cfa[indx - w1] - 18.0f * cfa[indx] * cfa[indx + w1] - 36.0f * cfa[indx] * cfa[indx - w2] - 36.0f * cfa[indx] * cfa[indx + w2] + 18.0f * cfa[indx] * cfa[indx - w3] + 18.0f * cfa[indx] * cfa[indx + w3] - 2.0f * cfa[indx] * cfa[indx - w4] - 2.0f * cfa[indx] * cfa[indx + w4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx - w1] * cfa[indx + w1] - 12.0f * cfa[indx - w1] * cfa[indx - w2] + 24.0f * cfa[indx - w1] * cfa[indx + w2] - 38.0f * cfa[indx - w1] * cfa[indx - w3] + 16.0f * cfa[indx - w1] * cfa[indx + w3] + 12.0f * cfa[indx - w1] * cfa[indx - w4] - 6.0f * cfa[indx - w1] * cfa[indx + w4] + 46.0f * cfa[indx - w1] * cfa[indx - w1] + 24.0f * cfa[indx + w1] * cfa[indx - w2] - 12.0f * cfa[indx + w1] * cfa[indx + w2] + 16.0f * cfa[indx + w1] * cfa[indx - w3] - 38.0f * cfa[indx + w1] * cfa[indx + w3] - 6.0f * cfa[indx + w1] * cfa[indx - w4] + 12.0f * cfa[indx + w1] * cfa[indx + w4] + 46.0f * cfa[indx + w1] * cfa[indx + w1] + 14.0f * cfa[indx - w2] * cfa[indx + w2] - 12.0f * cfa[indx - w2] * cfa[indx + w3] - 2.0f * cfa[indx - w2] * cfa[indx - w4] + 2.0f * cfa[indx - w2] * cfa[indx + w4] + 11.f * cfa[indx - w2] * cfa[indx - w2] - 12.0f * cfa[indx + w2] * cfa[indx - w3] + 2.0f * cfa[indx + w2] * cfa[indx - w4] - 2.0f * cfa[indx + w2] * cfa[indx + w4] + 11.f * cfa[indx + w2] * cfa[indx + w2] + 2.0f * cfa[indx - w3] * cfa[indx + w3] - 6.0f * cfa[indx - w3] * cfa[indx - w4] + 10.0f * cfa[indx - w3] * cfa[indx - w3] - 6.0f * cfa[indx + w3] * cfa[indx + w4] + 10.0f * cfa[indx + w3] * cfa[indx + w3] + 1.f * cfa[indx - w4] * cfa[indx - w4] + 1.f * cfa[indx + w4] * cfa[indx + w4];
            H_Stat = epssq - 18.0f  *  cfa[indx] * cfa[indx -  1] - 18.0f * cfa[indx] * cfa[indx +  1] - 36.0f * cfa[indx] * cfa[indx -  2] - 36.0f * cfa[indx] * cfa[indx +  2] + 18.0f * cfa[indx] * cfa[indx -  3] + 18.0f * cfa[indx] * cfa[indx +  3] - 2.0f * cfa[indx] * cfa[indx -  4] - 2.0f * cfa[indx] * cfa[indx +  4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx -  1] * cfa[indx +  1] - 12.0f * cfa[indx -  1] * cfa[indx -  2] + 24.0f * cfa[indx -  1] * cfa[indx +  2] - 38.0f * cfa[indx -  1] * cfa[indx -  3] + 16.0f * cfa[indx -  1] * cfa[indx +  3] + 12.0f * cfa[indx -  1] * cfa[indx -  4] - 6.0f * cfa[indx -  1] * cfa[indx +  4] + 46.0f * cfa[indx -  1] * cfa[indx -  1] + 24.0f * cfa[indx +  1] * cfa[indx -  2] - 12.0f * cfa[indx +  1] * cfa[indx +  2] + 16.0f * cfa[indx +  1] * cfa[indx -  3] - 38.0f * cfa[indx +  1] * cfa[indx +  3] - 6.0f * cfa[indx +  1] * cfa[indx -  4] + 12.0f * cfa[indx +  1] * cfa[indx +  4] + 46.0f * cfa[indx +  1] * cfa[indx +  1] + 14.0f * cfa[indx -  2] * cfa[indx +  2] - 12.0f * cfa[indx -  2] * cfa[indx +  3] - 2.0f * cfa[indx -  2] * cfa[indx -  4] + 2.0f * cfa[indx -  2] * cfa[indx +  4] + 11.f * cfa[indx -  2] * cfa[indx -  2] - 12.0f * cfa[indx +  2] * cfa[indx -  3] + 2.0f * cfa[indx +  2] * cfa[indx -  4] - 2.0f * cfa[indx +  2] * cfa[indx +  4] + 11.f * cfa[indx +  2] * cfa[indx +  2] + 2.0f * cfa[indx -  3] * cfa[indx +  3] - 6.0f * cfa[indx -  3] * cfa[indx -  4] + 10.0f * cfa[indx -  3] * cfa[indx -  3] - 6.0f * cfa[indx +  3] * cfa[indx +  4] + 10.0f * cfa[indx +  3] * cfa[indx +  3] + 1.f * cfa[indx -  4] * cfa[indx -  4] + 1.f * cfa[indx +  4] * cfa[indx +  4];

            VH_Dir[indx] = V_Stat / (V_Stat + H_Stat);

            //Calculate P/Q local discrimination
            P_Stat = epssq - 18.0f * cfa[indx] * cfa[indx - w1 - 1] - 18.0f * cfa[indx] * cfa[indx + w1 + 1] - 36.0f * cfa[indx] * cfa[indx - w2 - 2] - 36.0f * cfa[indx] * cfa[indx + w2 + 2] + 18.0f * cfa[indx] * cfa[indx - w3 - 3] + 18.0f * cfa[indx] * cfa[indx + w3 + 3] - 2.0f * cfa[indx] * cfa[indx - w4 - 4] - 2.0f * cfa[indx] * cfa[indx + w4 + 4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx - w1 - 1] * cfa[indx + w1 + 1] - 12.0f * cfa[indx - w1 - 1] * cfa[indx - w2 - 2] + 24.0f * cfa[indx - w1 - 1] * cfa[indx + w2 + 2] - 38.0f * cfa[indx - w1 - 1] * cfa[indx - w3 - 3] + 16.0f * cfa[indx - w1 - 1] * cfa[indx + w3 + 3] + 12.0f * cfa[indx - w1 - 1] * cfa[indx - w4 - 4] - 6.0f * cfa[indx - w1 - 1] * cfa[indx + w4 + 4] + 46.0f * cfa[indx - w1 - 1] * cfa[indx - w1 - 1] + 24.0f * cfa[indx + w1 + 1] * cfa[indx - w2 - 2] - 12.0f * cfa[indx + w1 + 1] * cfa[indx + w2 + 2] + 16.0f * cfa[indx + w1 + 1] * cfa[indx - w3 - 3] - 38.0f * cfa[indx + w1 + 1] * cfa[indx + w3 + 3] - 6.0f * cfa[indx + w1 + 1] * cfa[indx - w4 - 4] + 12.0f * cfa[indx + w1 + 1] * cfa[indx + w4 + 4] + 46.0f * cfa[indx + w1 + 1] * cfa[indx + w1 + 1] + 14.0f * cfa[indx - w2 - 2] * cfa[indx + w2 + 2] - 12.0f * cfa[indx - w2 - 2] * cfa[indx + w3 + 3] - 2.0f * cfa[indx - w2 - 2] * cfa[indx - w4 - 4] + 2.0f * cfa[indx - w2 - 2] * cfa[indx + w4 + 4] + 11.f * cfa[indx - w2 - 2] * cfa[indx - w2 - 2] - 12.0f * cfa[indx + w2 + 2] * cfa[indx - w3 - 3] + 2 * cfa[indx + w2 + 2] * cfa[indx - w4 - 4] - 2.0f * cfa[indx + w2 + 2] * cfa[indx + w4 + 4] + 11.f * cfa[indx + w2 + 2] * cfa[indx + w2 + 2] + 2.0f * cfa[indx - w3 - 3] * cfa[indx + w3 + 3] - 6.0f * cfa[indx - w3 - 3] * cfa[indx - w4 - 4] + 10.0f * cfa[indx - w3 - 3] * cfa[indx - w3 - 3] - 6.0f * cfa[indx + w3 + 3] * cfa[indx + w4 + 4] + 10.0f * cfa[indx + w3 + 3] * cfa[indx + w3 + 3] + 1.f * cfa[indx - w4 - 4] * cfa[indx - w4 - 4] + 1.f * cfa[indx + w4 + 4] * cfa[indx + w4 + 4];
            Q_Stat = epssq - 18.0f * cfa[indx] * cfa[indx + w1 - 1] - 18.0f * cfa[indx] * cfa[indx - w1 + 1] - 36.0f * cfa[indx] * cfa[indx + w2 - 2] - 36.0f * cfa[indx] * cfa[indx - w2 + 2] + 18.0f * cfa[indx] * cfa[indx + w3 - 3] + 18.0f * cfa[indx] * cfa[indx - w3 + 3] - 2.0f * cfa[indx] * cfa[indx + w4 - 4] - 2.0f * cfa[indx] * cfa[indx - w4 + 4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx + w1 - 1] * cfa[indx - w1 + 1] - 12.0f * cfa[indx + w1 - 1] * cfa[indx + w2 - 2] + 24.0f * cfa[indx + w1 - 1] * cfa[indx - w2 + 2] - 38.0f * cfa[indx + w1 - 1] * cfa[indx + w3 - 3] + 16.0f * cfa[indx + w1 - 1] * cfa[indx - w3 + 3] + 12.0f * cfa[indx + w1 - 1] * cfa[indx + w4 - 4] - 6.0f * cfa[indx + w1 - 1] * cfa[indx - w4 + 4] + 46.0f * cfa[indx + w1 - 1] * cfa[indx + w1 - 1] + 24.0f * cfa[indx - w1 + 1] * cfa[indx + w2 - 2] - 12.0f * cfa[indx - w1 + 1] * cfa[indx - w2 + 2] + 16.0f * cfa[indx - w1 + 1] * cfa[indx + w3 - 3] - 38.0f * cfa[indx - w1 + 1] * cfa[indx - w3 + 3] - 6.0f * cfa[indx - w1 + 1] * cfa[indx + w4 - 4] + 12.0f * cfa[indx - w1 + 1] * cfa[indx - w4 + 4] + 46.0f * cfa[indx - w1 + 1] * cfa[indx - w1 + 1] + 14.0f * cfa[indx + w2 - 2] * cfa[indx - w2 + 2] - 12.0f * cfa[indx + w2 - 2] * cfa[indx - w3 + 3] - 2.0f * cfa[indx + w2 - 2] * cfa[indx + w4 - 4] + 2.0f * cfa[indx + w2 - 2] * cfa[indx - w4 + 4] + 11.f * cfa[indx + w2 - 2] * cfa[indx + w2 - 2] - 12.0f * cfa[indx - w2 + 2] * cfa[indx + w3 - 3] + 2 * cfa[indx - w2 + 2] * cfa[indx + w4 - 4] - 2.0f * cfa[indx - w2 + 2] * cfa[indx - w4 + 4] + 11.f * cfa[indx - w2 + 2] * cfa[indx - w2 + 2] + 2.0f * cfa[indx + w3 - 3] * cfa[indx - w3 + 3] - 6.0f * cfa[indx + w3 - 3] * cfa[indx + w4 - 4] + 10.0f * cfa[indx + w3 - 3] * cfa[indx + w3 - 3] - 6.0f * cfa[indx - w3 + 3] * cfa[indx - w4 + 4] + 10.0f * cfa[indx - w3 + 3] * cfa[indx - w3 + 3] + 1.f * cfa[indx + w4 - 4] * cfa[indx + w4 - 4] + 1.f * cfa[indx - w4 + 4] * cfa[indx - w4 + 4];

            PQ_Dir[indx] = P_Stat / ( P_Stat + Q_Stat );

        }
    }

    VH_Disc = ( float ( * ) ) calloc( width * height, sizeof *VH_Disc ); merror ( VH_Disc, "rcd_demosaicing_171117()" );
    PQ_Disc = ( float ( * ) ) calloc( width * height, sizeof *PQ_Disc ); merror ( PQ_Disc, "rcd_demosaicing_171117()" );

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4, indx = row * width + col; col < width - 4; col++, indx++ ) {

            //Refined V/H local discrimination
            VH_Central_Value   = VH_Dir[indx];
            VH_Neighbour_Value = 0.25f * (VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1] + VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1]);

            VH_Disc[indx] = ( fabs( 0.5f - VH_Central_Value ) < fabs( 0.5f - VH_Neighbour_Value ) ) ? VH_Neighbour_Value : VH_Central_Value;

            //Refined P/Q local discrimination
            PQ_Central_Value   = PQ_Dir[indx];
            PQ_Neighbour_Value = 0.25f * (PQ_Dir[indx - w1 - 1] + PQ_Dir[indx - w1 + 1] + PQ_Dir[indx + w1 - 1] + PQ_Dir[indx + w1 + 1]);

            PQ_Disc[indx] = ( fabs( 0.5f - PQ_Central_Value ) < fabs( 0.5f - PQ_Neighbour_Value ) ) ? PQ_Neighbour_Value : PQ_Central_Value;

        }
    }

    free( VH_Dir );
    free( PQ_Dir );


    /**
    * STEP 2: Calculate the low pass filter
    */

    lpf = ( float ( * ) ) calloc( width * height, sizeof *lpf ); merror ( lpf, "rcd_demosaicing_171117()" );

    for ( row = 2; row < height - 2; row++ ) {
        for ( col = 2, indx = row * width + col; col < width - 2; col++, indx++ ) {

            //Low pass filter incorporating red and blue local samples
            lpf[indx] = 0.25f * cfa[indx] + 0.125f * ( cfa[indx - w1] + cfa[indx + w1] + cfa[indx - 1] + cfa[indx + 1] ) + 0.0625f * ( cfa[indx - w1 - 1] + cfa[indx - w1 + 1] + cfa[indx + w1 - 1] + cfa[indx + w1 + 1] );

        }
    }


    /**
    * STEP 3: Populate the green channel
    */

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + ( FC( row, 0 )&1 ), indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            //Cardinal gradients
            N_Grad = eps + fabs( cfa[indx - w1] - cfa[indx + w1] ) + fabs( cfa[indx] - cfa[indx - w2] ) + fabs( cfa[indx - w1] - cfa[indx - w3] ) + fabs( cfa[indx - w2] - cfa[indx - w4] );
            S_Grad = eps + fabs( cfa[indx + w1] - cfa[indx - w1] ) + fabs( cfa[indx] - cfa[indx + w2] ) + fabs( cfa[indx + w1] - cfa[indx + w3] ) + fabs( cfa[indx + w2] - cfa[indx + w4] );
            W_Grad = eps + fabs( cfa[indx -  1] - cfa[indx +  1] ) + fabs( cfa[indx] - cfa[indx -  2] ) + fabs( cfa[indx -  1] - cfa[indx -  3] ) + fabs( cfa[indx -  2] - cfa[indx -  4] );
            E_Grad = eps + fabs( cfa[indx +  1] - cfa[indx -  1] ) + fabs( cfa[indx] - cfa[indx +  2] ) + fabs( cfa[indx +  1] - cfa[indx +  3] ) + fabs( cfa[indx +  2] - cfa[indx +  4] );

            //Cardinal pixel estimations
            N_Est = cfa[indx - w1] * ( 1.f + ( lpf[indx] - lpf[indx - w2] ) / ( eps + lpf[indx] + lpf[indx - w2] ) );
            S_Est = cfa[indx + w1] * ( 1.f + ( lpf[indx] - lpf[indx + w2] ) / ( eps + lpf[indx] + lpf[indx + w2] ) );
            W_Est = cfa[indx -  1] * ( 1.f + ( lpf[indx] - lpf[indx -  2] ) / ( eps + lpf[indx] + lpf[indx -  2] ) );
            E_Est = cfa[indx +  1] * ( 1.f + ( lpf[indx] - lpf[indx +  2] ) / ( eps + lpf[indx] + lpf[indx +  2] ) );

            //Interpolate G@R & G@B
            V_Est = ( S_Grad * N_Est + N_Grad * S_Est ) / ( N_Grad + S_Grad );
            H_Est = ( W_Grad * E_Est + E_Grad * W_Est ) / ( E_Grad + W_Grad );

            rgb[indx][1] = LIM( VH_Disc[indx] * H_Est + ( 1.f - VH_Disc[indx] ) * V_Est, 0.f, 1.f );

        }
    }

    free( lpf );
    free( cfa );


    /**
    * STEP 4: Populate the red and blue channel
    */

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + ( FC( row, 0 )&1 ), indx = row * width + col, c = 2 - FC( row, col ); col < width - 4; col += 2, indx += 2 ) {

            //Diagonal gradients
            NW_Grad = eps + fabs( rgb[indx - w1 - 1][c] - rgb[indx + w1 + 1][c] ) + fabs( rgb[indx - w1 - 1][c] - rgb[indx - w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 - 2][1] );
            NE_Grad = eps + fabs( rgb[indx - w1 + 1][c] - rgb[indx + w1 - 1][c] ) + fabs( rgb[indx - w1 + 1][c] - rgb[indx - w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 + 2][1] );
            SW_Grad = eps + fabs( rgb[indx + w1 - 1][c] - rgb[indx - w1 + 1][c] ) + fabs( rgb[indx + w1 - 1][c] - rgb[indx + w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 - 2][1] );
            SE_Grad = eps + fabs( rgb[indx + w1 + 1][c] - rgb[indx - w1 - 1][c] ) + fabs( rgb[indx + w1 + 1][c] - rgb[indx + w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 + 2][1] );

            //Diagonal colour differences
            NW_Est = rgb[indx - w1 - 1][c] - rgb[indx - w1 - 1][1];
            NE_Est = rgb[indx - w1 + 1][c] - rgb[indx - w1 + 1][1];
            SW_Est = rgb[indx + w1 - 1][c] - rgb[indx + w1 - 1][1];
            SE_Est = rgb[indx + w1 + 1][c] - rgb[indx + w1 + 1][1];

            //Interpolate R@B and B@R
            P_Est = ( NW_Grad * SE_Est + SE_Grad * NW_Est ) / ( NW_Grad + SE_Grad );
            Q_Est = ( NE_Grad * SW_Est + SW_Grad * NE_Est ) / ( NE_Grad + SW_Grad );

            rgb[indx][c] = LIM( rgb[indx][1] + ( 1.f - PQ_Disc[indx] ) * P_Est + PQ_Disc[indx] * Q_Est, 0.f, 1.f );

        }
    }

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + ( FC( row, 1 )&1 ), indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            for ( c = 0; c <= 2; c += 2 ) {

                //Cardinal gradients
                N_Grad = eps + fabs( rgb[indx][1] - rgb[indx - w2][1] ) + fabs( rgb[indx - w1][c] - rgb[indx + w1][c] ) + fabs( rgb[indx - w1][c] - rgb[indx - w3][c] );
                S_Grad = eps + fabs( rgb[indx][1] - rgb[indx + w2][1] ) + fabs( rgb[indx + w1][c] - rgb[indx - w1][c] ) + fabs( rgb[indx + w1][c] - rgb[indx + w3][c] );
                W_Grad = eps + fabs( rgb[indx][1] - rgb[indx -  2][1] ) + fabs( rgb[indx -  1][c] - rgb[indx +  1][c] ) + fabs( rgb[indx -  1][c] - rgb[indx -  3][c] );
                E_Grad = eps + fabs( rgb[indx][1] - rgb[indx +  2][1] ) + fabs( rgb[indx +  1][c] - rgb[indx -  1][c] ) + fabs( rgb[indx +  1][c] - rgb[indx +  3][c] );

                //Cardinal colour differences
                N_Est = rgb[indx - w1][c] - rgb[indx - w1][1];
                S_Est = rgb[indx + w1][c] - rgb[indx + w1][1];
                W_Est = rgb[indx -  1][c] - rgb[indx -  1][1];
                E_Est = rgb[indx +  1][c] - rgb[indx +  1][1];

                //Interpolate R@G and B@G
                V_Est = ( N_Grad * S_Est + S_Grad * N_Est ) / ( N_Grad + S_Grad );
                H_Est = ( E_Grad * W_Est + W_Grad * E_Est ) / ( E_Grad + W_Grad );

                rgb[indx][c] = LIM( rgb[indx][1] + ( 1.f - VH_Disc[indx] ) * V_Est + VH_Disc[indx] * H_Est, 0.f, 1.f );

            }
        }
    }

    free( PQ_Disc );
    free( VH_Disc );


    //Revert back the floating point buffers to dcraw *image
    for ( row = 0; row < height; row++ ) {
        for ( col = 0, indx = row * width + col; col < width; col++, indx++ ) {

            FORC3 image[indx][c] = (ushort)CLIP( 65535.f * rgb[indx][c] );

        }
    }

    border_interpolate(4);

    free( rgb );

}