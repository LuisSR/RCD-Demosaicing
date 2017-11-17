/**
* RATIO CORRECTED DEMOSAICING
* Luis Sanz Rodriguez (luis.sanz.rodriguez(at)gmail(dot)com)
*
* Release 2.2 @ 171117
*/

void CLASS rcd_demosaicing_171117()
{
    int row, col, indx, c;
    int w1 = width, w2 = 2 * width, w3 = 3 * width, w4 = 4 * width;

    //Tolerance to avoid dividing by zero
    static const float eps = 1e-5, epssq = 1e-10;

    //Gradients
    float Ngrad, Egrad, Wgrad, Sgrad, NWgrad, NEgrad, SWgrad, SEgrad;

    //Pixel estimation
    float Nest, Eest, West, Sest, NWest, NEest, SWest, SEest;
    float vint, hint, pint, qint;

    //Directional discrimination
    float vstat, hstat, vstatalt, mstat, pstat, mstatalt;
    float vh_central_value, vh_neighbour_value, mp_central_value, mp_neighbour_value;
    float ( *vhdir ), ( *vhdisc ), ( *mpdir ), ( *mpdisc );

    //Low pass filter
    float ( *lpf );


    /**
    * STEP 1: Find cardinal and diagonal interpolation directions
    */

    vhdir = ( float ( * ) ) calloc( width * height, sizeof *vhdir ); merror ( vhdir, "rcd_demosaicing_171117()" );
    mpdir = ( float ( * ) ) calloc( width * height, sizeof *mpdir ); merror ( mpdir, "rcd_demosaicing_171117()" );

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4, indx = row * width + col; col < width - 4; col++, indx++ ) {

            //Calculate h/v local discrimination
            vstat = epssq - 18.0f  *  cfa[indx] * cfa[indx - w1] - 18.0f * cfa[indx] * cfa[indx + w1] - 36.0f * cfa[indx] * cfa[indx - w2] - 36.0f * cfa[indx] * cfa[indx + w2] + 18.0f * cfa[indx] * cfa[indx - w3] + 18.0f * cfa[indx] * cfa[indx + w3] - 2.0f * cfa[indx] * cfa[indx - w4] - 2.0f * cfa[indx] * cfa[indx + w4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx - w1] * cfa[indx + w1] - 12.0f * cfa[indx - w1] * cfa[indx - w2] + 24.0f * cfa[indx - w1] * cfa[indx + w2] - 38.0f * cfa[indx - w1] * cfa[indx - w3] + 16.0f * cfa[indx - w1] * cfa[indx + w3] + 12.0f * cfa[indx - w1] * cfa[indx - w4] - 6.0f * cfa[indx - w1] * cfa[indx + w4] + 46.0f * cfa[indx - w1] * cfa[indx - w1] + 24.0f * cfa[indx + w1] * cfa[indx - w2] - 12.0f * cfa[indx + w1] * cfa[indx + w2] + 16.0f * cfa[indx + w1] * cfa[indx - w3] - 38.0f * cfa[indx + w1] * cfa[indx + w3] - 6.0f * cfa[indx + w1] * cfa[indx - w4] + 12.0f * cfa[indx + w1] * cfa[indx + w4] + 46.0f * cfa[indx + w1] * cfa[indx + w1] + 14.0f * cfa[indx - w2] * cfa[indx + w2] - 12.0f * cfa[indx - w2] * cfa[indx + w3] - 2.0f * cfa[indx - w2] * cfa[indx - w4] + 2.0f * cfa[indx - w2] * cfa[indx + w4] + 11.0f * cfa[indx - w2] * cfa[indx - w2] - 12.0f * cfa[indx + w2] * cfa[indx - w3] + 2.0f * cfa[indx + w2] * cfa[indx - w4] - 2.0f * cfa[indx + w2] * cfa[indx + w4] + 11.0f * cfa[indx + w2] * cfa[indx + w2] + 2.0f * cfa[indx - w3] * cfa[indx + w3] - 6.0f * cfa[indx - w3] * cfa[indx - w4] + 10.0f * cfa[indx - w3] * cfa[indx - w3] - 6.0f * cfa[indx + w3] * cfa[indx + w4] + 10.0f * cfa[indx + w3] * cfa[indx + w3] + 1.0f * cfa[indx - w4] * cfa[indx - w4] + 1.0f * cfa[indx + w4] * cfa[indx + w4];
            hstat = epssq - 18.0f  *  cfa[indx] * cfa[indx -  1] - 18.0f * cfa[indx] * cfa[indx +  1] - 36.0f * cfa[indx] * cfa[indx -  2] - 36.0f * cfa[indx] * cfa[indx +  2] + 18.0f * cfa[indx] * cfa[indx -  3] + 18.0f * cfa[indx] * cfa[indx +  3] - 2.0f * cfa[indx] * cfa[indx -  4] - 2.0f * cfa[indx] * cfa[indx +  4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx -  1] * cfa[indx +  1] - 12.0f * cfa[indx -  1] * cfa[indx -  2] + 24.0f * cfa[indx -  1] * cfa[indx +  2] - 38.0f * cfa[indx -  1] * cfa[indx -  3] + 16.0f * cfa[indx -  1] * cfa[indx +  3] + 12.0f * cfa[indx -  1] * cfa[indx -  4] - 6.0f * cfa[indx -  1] * cfa[indx +  4] + 46.0f * cfa[indx -  1] * cfa[indx -  1] + 24.0f * cfa[indx +  1] * cfa[indx -  2] - 12.0f * cfa[indx +  1] * cfa[indx +  2] + 16.0f * cfa[indx +  1] * cfa[indx -  3] - 38.0f * cfa[indx +  1] * cfa[indx +  3] - 6.0f * cfa[indx +  1] * cfa[indx -  4] + 12.0f * cfa[indx +  1] * cfa[indx +  4] + 46.0f * cfa[indx +  1] * cfa[indx +  1] + 14.0f * cfa[indx -  2] * cfa[indx +  2] - 12.0f * cfa[indx -  2] * cfa[indx +  3] - 2.0f * cfa[indx -  2] * cfa[indx -  4] + 2.0f * cfa[indx -  2] * cfa[indx +  4] + 11.0f * cfa[indx -  2] * cfa[indx -  2] - 12.0f * cfa[indx +  2] * cfa[indx -  3] + 2.0f * cfa[indx +  2] * cfa[indx -  4] - 2.0f * cfa[indx +  2] * cfa[indx +  4] + 11.0f * cfa[indx +  2] * cfa[indx +  2] + 2.0f * cfa[indx -  3] * cfa[indx +  3] - 6.0f * cfa[indx -  3] * cfa[indx -  4] + 10.0f * cfa[indx -  3] * cfa[indx -  3] - 6.0f * cfa[indx +  3] * cfa[indx +  4] + 10.0f * cfa[indx +  3] * cfa[indx +  3] + 1.0f * cfa[indx -  4] * cfa[indx -  4] + 1.0f * cfa[indx +  4] * cfa[indx +  4];

            vhdir[indx] = vstat / (vstat + hstat);

            //Calculate m/p local discrimination
            mstat = epssq - 18.0f * cfa[indx] * cfa[indx - w1 - 1] - 18.0f * cfa[indx] * cfa[indx + w1 + 1] - 36.0f * cfa[indx] * cfa[indx - w2 - 2] - 36.0f * cfa[indx] * cfa[indx + w2 + 2] + 18.0f * cfa[indx] * cfa[indx - w3 - 3] + 18.0f * cfa[indx] * cfa[indx + w3 + 3] - 2.0f * cfa[indx] * cfa[indx - w4 - 4] - 2.0f * cfa[indx] * cfa[indx + w4 + 4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx - w1 - 1] * cfa[indx + w1 + 1] - 12.0f * cfa[indx - w1 - 1] * cfa[indx - w2 - 2] + 24.0f * cfa[indx - w1 - 1] * cfa[indx + w2 + 2] - 38.0f * cfa[indx - w1 - 1] * cfa[indx - w3 - 3] + 16.0f * cfa[indx - w1 - 1] * cfa[indx + w3 + 3] + 12.0f * cfa[indx - w1 - 1] * cfa[indx - w4 - 4] - 6.0f * cfa[indx - w1 - 1] * cfa[indx + w4 + 4] + 46.0f * cfa[indx - w1 - 1] * cfa[indx - w1 - 1] + 24.0f * cfa[indx + w1 + 1] * cfa[indx - w2 - 2] - 12.0f * cfa[indx + w1 + 1] * cfa[indx + w2 + 2] + 16.0f * cfa[indx + w1 + 1] * cfa[indx - w3 - 3] - 38.0f * cfa[indx + w1 + 1] * cfa[indx + w3 + 3] - 6.0f * cfa[indx + w1 + 1] * cfa[indx - w4 - 4] + 12.0f * cfa[indx + w1 + 1] * cfa[indx + w4 + 4] + 46.0f * cfa[indx + w1 + 1] * cfa[indx + w1 + 1] + 14.0f * cfa[indx - w2 - 2] * cfa[indx + w2 + 2] - 12.0f * cfa[indx - w2 - 2] * cfa[indx + w3 + 3] - 2.0f * cfa[indx - w2 - 2] * cfa[indx - w4 - 4] + 2.0f * cfa[indx - w2 - 2] * cfa[indx + w4 + 4] + 11.0f * cfa[indx - w2 - 2] * cfa[indx - w2 - 2] - 12.0f * cfa[indx + w2 + 2] * cfa[indx - w3 - 3] + 2 * cfa[indx + w2 + 2] * cfa[indx - w4 - 4] - 2.0f * cfa[indx + w2 + 2] * cfa[indx + w4 + 4] + 11.0f * cfa[indx + w2 + 2] * cfa[indx + w2 + 2] + 2.0f * cfa[indx - w3 - 3] * cfa[indx + w3 + 3] - 6.0f * cfa[indx - w3 - 3] * cfa[indx - w4 - 4] + 10.0f * cfa[indx - w3 - 3] * cfa[indx - w3 - 3] - 6.0f * cfa[indx + w3 + 3] * cfa[indx + w4 + 4] + 10.0f * cfa[indx + w3 + 3] * cfa[indx + w3 + 3] + 1.0f * cfa[indx - w4 - 4] * cfa[indx - w4 - 4] + 1.0f * cfa[indx + w4 + 4] * cfa[indx + w4 + 4];
            pstat = epssq - 18.0f * cfa[indx] * cfa[indx + w1 - 1] - 18.0f * cfa[indx] * cfa[indx - w1 + 1] - 36.0f * cfa[indx] * cfa[indx + w2 - 2] - 36.0f * cfa[indx] * cfa[indx - w2 + 2] + 18.0f * cfa[indx] * cfa[indx + w3 - 3] + 18.0f * cfa[indx] * cfa[indx - w3 + 3] - 2.0f * cfa[indx] * cfa[indx + w4 - 4] - 2.0f * cfa[indx] * cfa[indx - w4 + 4] + 38.0f * cfa[indx] * cfa[indx] - 70.0f * cfa[indx + w1 - 1] * cfa[indx - w1 + 1] - 12.0f * cfa[indx + w1 - 1] * cfa[indx + w2 - 2] + 24.0f * cfa[indx + w1 - 1] * cfa[indx - w2 + 2] - 38.0f * cfa[indx + w1 - 1] * cfa[indx + w3 - 3] + 16.0f * cfa[indx + w1 - 1] * cfa[indx - w3 + 3] + 12.0f * cfa[indx + w1 - 1] * cfa[indx + w4 - 4] - 6.0f * cfa[indx + w1 - 1] * cfa[indx - w4 + 4] + 46.0f * cfa[indx + w1 - 1] * cfa[indx + w1 - 1] + 24.0f * cfa[indx - w1 + 1] * cfa[indx + w2 - 2] - 12.0f * cfa[indx - w1 + 1] * cfa[indx - w2 + 2] + 16.0f * cfa[indx - w1 + 1] * cfa[indx + w3 - 3] - 38.0f * cfa[indx - w1 + 1] * cfa[indx - w3 + 3] - 6.0f * cfa[indx - w1 + 1] * cfa[indx + w4 - 4] + 12.0f * cfa[indx - w1 + 1] * cfa[indx - w4 + 4] + 46.0f * cfa[indx - w1 + 1] * cfa[indx - w1 + 1] + 14.0f * cfa[indx + w2 - 2] * cfa[indx - w2 + 2] - 12.0f * cfa[indx + w2 - 2] * cfa[indx - w3 + 3] - 2.0f * cfa[indx + w2 - 2] * cfa[indx + w4 - 4] + 2.0f * cfa[indx + w2 - 2] * cfa[indx - w4 + 4] + 11.0f * cfa[indx + w2 - 2] * cfa[indx + w2 - 2] - 12.0f * cfa[indx - w2 + 2] * cfa[indx + w3 - 3] + 2 * cfa[indx - w2 + 2] * cfa[indx + w4 - 4] - 2.0f * cfa[indx - w2 + 2] * cfa[indx - w4 + 4] + 11.0f * cfa[indx - w2 + 2] * cfa[indx - w2 + 2] + 2.0f * cfa[indx + w3 - 3] * cfa[indx - w3 + 3] - 6.0f * cfa[indx + w3 - 3] * cfa[indx + w4 - 4] + 10.0f * cfa[indx + w3 - 3] * cfa[indx + w3 - 3] - 6.0f * cfa[indx - w3 + 3] * cfa[indx - w4 + 4] + 10.0f * cfa[indx - w3 + 3] * cfa[indx - w3 + 3] + 1.0f * cfa[indx + w4 - 4] * cfa[indx + w4 - 4] + 1.0f * cfa[indx - w4 + 4] * cfa[indx - w4 + 4];

            mpdir[indx] = mstat / ( mstat + pstat );

        }
    }

    vhdisc = ( float ( * ) ) calloc( width * height, sizeof *vhdisc ); merror ( vhdisc, "rcd_demosaicing_171117()" );
    mpdisc = ( float ( * ) ) calloc( width * height, sizeof *mpdisc ); merror ( mpdisc, "rcd_demosaicing_171117()" );

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4, indx = row * width + col; col < width - 4; col++, indx++ ) {

            //Refined h/v local discrimination
            vh_central_value   = vhdir[indx];
            vh_neighbour_value = 0.25f * (vhdir[indx - w1 - 1] + vhdir[indx - w1 + 1] + vhdir[indx + w1 - 1] + vhdir[indx + w1 + 1]);

            vhdisc[indx] = ( fabs( 0.5f - vh_central_value ) < fabs( 0.5f - vh_neighbour_value ) ) ? vh_neighbour_value : vh_central_value;

            //Refined m/p local discrimination
            mp_central_value   = mpdir[indx];
            mp_neighbour_value = 0.25f * (mpdir[indx - w1 - 1] + mpdir[indx - w1 + 1] + mpdir[indx + w1 - 1] + mpdir[indx + w1 + 1]);

            mpdisc[indx] = ( fabs( 0.5f - mp_central_value ) < fabs( 0.5f - mp_neighbour_value ) ) ? mp_neighbour_value : mp_central_value;

        }
    }

    free( vhdir );
    free( mpdir );


    /**
    * STEP 2: Calculate the low pass filter
    */

    lpf = ( float ( * ) ) calloc( width * height, sizeof *lpf ); merror ( lpf, "rcd_demosaicing_171117()" );

    for ( row = 1; row < height - 1; row++ ) {
        for ( col = 1, indx = row * width + col; col < width - 1; col++, indx++ ) {

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
            Ngrad = eps + fabs( cfa[indx - w1] - cfa[indx + w1] ) + fabs( cfa[indx] - cfa[indx - w2] ) + fabs( cfa[indx - w1] - cfa[indx - w3] ) + fabs( cfa[indx - w2] - cfa[indx - w4] );
            Sgrad = eps + fabs( cfa[indx + w1] - cfa[indx - w1] ) + fabs( cfa[indx] - cfa[indx + w2] ) + fabs( cfa[indx + w1] - cfa[indx + w3] ) + fabs( cfa[indx + w2] - cfa[indx + w4] );
            Wgrad = eps + fabs( cfa[indx -  1] - cfa[indx +  1] ) + fabs( cfa[indx] - cfa[indx -  2] ) + fabs( cfa[indx -  1] - cfa[indx -  3] ) + fabs( cfa[indx -  2] - cfa[indx -  4] );
            Egrad = eps + fabs( cfa[indx +  1] - cfa[indx -  1] ) + fabs( cfa[indx] - cfa[indx +  2] ) + fabs( cfa[indx +  1] - cfa[indx +  3] ) + fabs( cfa[indx +  2] - cfa[indx +  4] );

            //Cardinal pixel estimations
            Nest = cfa[indx - w1] * ( 1.f + ( lpf[indx] - lpf[indx - w2] ) / ( eps + lpf[indx] + lpf[indx - w2] ) );
            Sest = cfa[indx + w1] * ( 1.f + ( lpf[indx] - lpf[indx + w2] ) / ( eps + lpf[indx] + lpf[indx + w2] ) );
            West = cfa[indx -  1] * ( 1.f + ( lpf[indx] - lpf[indx -  2] ) / ( eps + lpf[indx] + lpf[indx -  2] ) );
            Eest = cfa[indx +  1] * ( 1.f + ( lpf[indx] - lpf[indx +  2] ) / ( eps + lpf[indx] + lpf[indx +  2] ) );

            //Interpolate G@R & G@B
            vint = ( Sgrad * Nest + Ngrad * Sest ) / ( Ngrad + Sgrad );
            hint = ( Wgrad * Eest + Egrad * West ) / ( Egrad + Wgrad );

            rgb[indx][1] = LIM( vhdisc[indx] * hint + ( 1.0f - vhdisc[indx] ) * vint, 0.f, 1.f );

        }
    }

    free( lpf );


    /**
    * STEP 4: Populate the red and blue channel
    */

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + ( FC( row, 0 )&1 ), indx = row * width + col, c = 2 - FC( row, col ); col < width - 4; col += 2, indx += 2 ) {

            //Diagonal gradients
            NWgrad = eps + fabs( rgb[indx - w1 - 1][c] - rgb[indx + w1 + 1][c] ) + fabs( rgb[indx - w1 - 1][c] - rgb[indx - w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 - 2][1] );
            NEgrad = eps + fabs( rgb[indx - w1 + 1][c] - rgb[indx + w1 - 1][c] ) + fabs( rgb[indx - w1 + 1][c] - rgb[indx - w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx - w2 + 2][1] );
            SWgrad = eps + fabs( rgb[indx + w1 - 1][c] - rgb[indx - w1 + 1][c] ) + fabs( rgb[indx + w1 - 1][c] - rgb[indx + w3 - 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 - 2][1] );
            SEgrad = eps + fabs( rgb[indx + w1 + 1][c] - rgb[indx - w1 - 1][c] ) + fabs( rgb[indx + w1 + 1][c] - rgb[indx + w3 + 3][c] ) + fabs( rgb[indx][1] - rgb[indx + w2 + 2][1] );

            //Diagonal colour differences
            NWest = rgb[indx - w1 - 1][c] - rgb[indx - w1 - 1][1];
            NEest = rgb[indx - w1 + 1][c] - rgb[indx - w1 + 1][1];
            SWest = rgb[indx + w1 - 1][c] - rgb[indx + w1 - 1][1];
            SEest = rgb[indx + w1 + 1][c] - rgb[indx + w1 + 1][1];

            //Interpolate R@B and B@R
            pint = ( NWgrad * SEest + SEgrad * NWest ) / ( NWgrad + SEgrad );
            qint = ( NEgrad * SWest + SWgrad * NEest ) / ( NEgrad + SWgrad );

            rgb[indx][c] = LIM( rgb[indx][1] + ( 1.0f - mpdisc[indx] ) * pint + mpdisc[indx] * qint, 0.f, 1.f );

        }
    }

    for ( row = 4; row < height - 4; row++ ) {
        for ( col = 4 + ( FC( row, 1 )&1 ), indx = row * width + col; col < width - 4; col += 2, indx += 2 ) {

            for ( c = 0; c <= 2; c += 2 ) {

                //Cardinal gradients
                Ngrad = eps + fabs( rgb[indx][1] - rgb[indx - w2][1] ) + fabs( rgb[indx - w1][c] - rgb[indx + w1][c] ) + fabs( rgb[indx - w1][c] - rgb[indx - w3][c] );
                Sgrad = eps + fabs( rgb[indx][1] - rgb[indx + w2][1] ) + fabs( rgb[indx + w1][c] - rgb[indx - w1][c] ) + fabs( rgb[indx + w1][c] - rgb[indx + w3][c] );
                Wgrad = eps + fabs( rgb[indx][1] - rgb[indx -  2][1] ) + fabs( rgb[indx -  1][c] - rgb[indx +  1][c] ) + fabs( rgb[indx -  1][c] - rgb[indx -  3][c] );
                Egrad = eps + fabs( rgb[indx][1] - rgb[indx +  2][1] ) + fabs( rgb[indx +  1][c] - rgb[indx -  1][c] ) + fabs( rgb[indx +  1][c] - rgb[indx +  3][c] );

                //Cardinal colour differences
                Nest = rgb[indx - w1][c] - rgb[indx - w1][1];
                Sest = rgb[indx + w1][c] - rgb[indx + w1][1];
                West = rgb[indx -  1][c] - rgb[indx -  1][1];
                Eest = rgb[indx +  1][c] - rgb[indx +  1][1];

                //Interpolate R@G and B@G
                vint = ( Ngrad * Sest + Sgrad * Nest ) / ( Ngrad + Sgrad );
                hint = ( Egrad * West + Wgrad * Eest ) / ( Egrad + Wgrad );

                rgb[indx][c] = LIM( rgb[indx][1] + ( 1.0f - vhdisc[indx] ) * vint + vhdisc[indx] * hint, 0.f, 1.f );

            }
        }
    }

    free( mpdisc );
    free( vhdisc );

}