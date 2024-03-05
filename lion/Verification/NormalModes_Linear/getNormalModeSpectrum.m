function [freqs, modes] = getNormalModeSpectrum(N)
%GETNORMALMODESPECTRA Calculates normal mode spectra for N ion chain
%  Calculates the normal mode spectrum for a chain of N ions. The table of
%  analytic results is taken from James, App Phys B 66 (2) 1998, where the
%  frequencies of motion in the z direction are listed. The frequency
%  returned is normalised with respect to the quantity \nu, which is
%  defined as the angular trap frequency.
%  
%  Syntax:
%   freqs = getNormalModeSpectrum(N)

switch N
    case 1
        eigenVal = 1;
        modes = 1;
    case 2
        eigenVal = [ 1 3 ];
        modes = [ ...
           +0.7071 +0.7071;
           -0.7071 +0.7071;
           ];
    case 3
        eigenVal = [ 1 3 5.8 ];
        modes = [ ...
           +0.5774 +0.5774 +0.5774;
           -0.7071 +0.0000 +0.7071;
           +0.4082 -0.8165 +0.4082;
           ];
    case 4
        eigenVal = [ 1 3 5.81 9.308 ];
        modes = [ ...
            +0.5000 +0.5000 +0.5000 +0.5000;
            -0.6742 -0.2132 +0.2132 +0.6742;
            +0.5000 -0.5000 -0.5000 +0.5000;
            -0.2132 +0.6742 -0.6742 +0.2132;
            ];
    case 5
        eigenVal = [ 1 3 5.818 9.332 13.47 ];
        modes = [ ...
            +0.4472 +0.4472 +0.4472 +0.4472 +0.4472;
            -0.6392 -0.3017 +0.0000 +0.3017 +0.6395;
            +0.5377 -0.2805 -0.5143 -0.2805 +0.5377;
            -0.3017 +0.6395 +0.0000 -0.6395 +0.3017;
            +0.1045 -0.4704 +0.7318 -0.4704 +0.1045;
            ];
    case 6
        eigenVal = [ 1 3 5.824 9.352 13.51 18.27 ];
        modes = [ ...
            +0.4082 +0.4082 +0.4082 +0.4082 +0.4082 +0.4082;
            -0.6080 -0.3433 -0.1118 +0.1118 +0.3433 +0.6080;
            -0.5531 +0.1332 +0.4199 +0.4199 +0.1332 -0.5531;
            +0.3577 -0.5431 -0.2778 +0.2778 +0.5431 -0.3577;
            +0.1655 -0.5618 +0.3963 +0.3963 -0.5618 +0.1655;
            -0.0490 +0.2954 -0.6406 +0.6406 -0.2954 +0.0490;
            ];            
    case 7
        eigenVal = [ 1 3 5.829 9.369 13.55 18.32 23.66 ];
        modes = [ ...
            +0.3780 +0.3780 +0.3780 +0.3780 +0.3780 +0.3780 +0.3780;
            -0.5801 -0.3636 -0.1768 +0.0000 +0.1768 +0.3636 +0.5801;
            -0.5579 +0.0310 +0.3213 +0.4111 +0.3213 +0.0310 -0.5579;
            -0.3952 +0.4450 +0.3818 +0.0000 -0.3818 -0.4450 +0.3952;
            -0.2130 +0.5714 -0.1199 -0.4769 -0.1199 +0.5714 -0.2130;
            +0.0851 -0.4121 +0.5683 +0.0000 -0.5683 +0.4121 -0.0851;
            +0.0222 -0.1723 +0.4894 -0.6787 +0.4894 -0.1723 +0.0222;
            ];            
    case 8
        eigenVal = [ 1 3 5.834 9.383 13.58 18.37 23.73 29.63 ];
    case 9
        eigenVal = [ 1 3 5.838 9.396 13.6 18.41 23.79 29.71 36.16 ];
    case 10
        eigenVal = [ 1 3 5.841 9.408 13.63 18.45 23.85 29.79 36.26 43.24 ];
    otherwise
        error('Unsupported value of N');
end

% see eq (23), which relates the frequency of mode p to eigenvalue p.
freqs = eigenVal.^0.5;

end

