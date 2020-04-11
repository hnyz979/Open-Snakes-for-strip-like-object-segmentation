function filtered=gaussian_filter(unfiltered,sigma_filt)
smask = fspecial('gaussian', ceil(3*sigma_filt), sigma_filt);
filtered = filter2(smask, unfiltered, 'same');