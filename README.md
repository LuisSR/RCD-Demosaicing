# RCD-Demosaicing
A Bayer CFA demosaicing routine aimed to reduce pixel overshooting.

## Brief overview

The Ratio Corrected Demosaicing algorithm intends to smooth the colour correction errors that are common in many other interpolation methods. To achieve this goal, the algorithm operates at two different levels:
* It uses a novel directional discrimination statistic which is invariant to chromatic aberrations.
* It estimates the pixel value applying a smoothing ratio to the spatial difference correction in a low-pass filter domain instead of in the more usual colour difference in the CFA domain.

The algorithm comprises four parts:

### 1. Local directional discrimination
The algorithm calculates the relative 1D local gradients directly on the Bayer data for the horizonal, vertical, diagonal P and diagonal Q directions.

It then refines the calculation trying to infere more discriminational strength from the neighbourhood wherever possible.

### 2. Low-pass filter
Taking advantadge of the structure of the Bayer CFA, a very simple and stable low-pass filter is used to smooth the colour intensity differences between the red and the blue channels.

    1 2 1
    2 4 2 * CFA
    1 2 1

This convolution results in a grey scale image which is softer than the input but is almost artifact free.

### 3. Green channel interpolation
Using the data from the steps #1 and #2, the green channel is updated using a pixel estimation based on ratios. Instead of the widespread Hamilton-Adams interpolation, the proposed ratio-corrected interpolation applies the colour correction taking into account the context of the intensity difference between spatially separated pixels. This way, in hard edges colour overshooting issues are under control.

Example for the East pixel value estimation:

    Hamilton-Adams  => GE_Est(0) = G(1) + ( R(0) - R(2) ) / 2
    Ratio-corrected => GE_Est(0) = G(1) * ( 1 + ( R(0) - R(2) ) / ( R(0) + R(2) ) )

To further smooth the output, the low-pass filter from step #2 is introduced in the formula:

    Ratio-corrected LPF => GE_Est(0) = G(1) * ( 1 + ( LPF(0) - LPF(2) ) / ( LPF(0) + LPF(2) ) )

### 4. Red and Blue channels interpolation
The routine used for interpolating the chrominance is fairly simple yet surprisingly effective. Using the full green channel from step #3, the red and blue channels are updated with the local colour differences.

The gradients are empirically optimized.

## Usage

I provide binary files built in Visual Studio C++ 2015 under Windows 8.1 (64-bit) and GCC under Lubuntu 16.04 (32-bit).

To run RCD, use the attribute "-q 4". Command example:

    dcraw -v -T -6 -H 0 -W -o 1 -w -q 4 file.raw

The algorithm can be compiled manually. In VSC++ 2015 I run

    cl /nologo /Ox /w /DWIN32 /DDJGPP /DNO_JASPER /DNO_LCMS /DNO_JPEG /Iinclude /Fedcraw.exe dcraw.c /link /subsystem:console user32.lib

In GCC, I use

    gcc -o dcraw -O4 ./dcraw.c -w -lm -DNO_JPEG -DNO_LCMS -DNO_JASPER
