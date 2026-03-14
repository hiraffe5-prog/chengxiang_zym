function[shift_row,shift_col]=pixel_registration(slc1,slc2)

    [rows,cols]=size(slc1);
    cross_corr = fftshift(ifft2(fft2(slc1,rows,cols).*conj(fft2(slc2,rows,cols))));%∆µ”Úª•œ‡πÿ
    [shift_row,shift_col]=find(abs(cross_corr)==max(abs(cross_corr(:))));
    shift_row = fix((rows+1)/2-shift_row);
    shift_col = fix((cols+1)/2-shift_col);