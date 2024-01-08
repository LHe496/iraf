include <mach.h>
include	<imhdr.h>
include <fset.h>
include	"vt.h"

define	WDSBRSTR	50

# READSS1 -- Read a type 1 sector scan from tape and format into 3 iraf images.
# Type one sector scans consist of three images packed into 32 bits per
# pixel.  The three images are 1. velocity (12 bits) 2. select (12 bits) and
# 3. continuum intensity (8 bits).  The images are only 256 pixels high as
# opposed to 512 pixels high for the other scans.

procedure readss1 (inputfd, filenumber, brief, select, bright, velocity, hs)

int	inputfd			# file descriptor for input (usually tape)
int	filenumber		# file number on tape
bool	brief			# short output file names
bool	select			# flag to make select image
bool	bright			# flag to make bright image
bool	velocity		# flag to make velocity image
int	hs			# header data structure pointer

char	velimage[SZ_FNAME]	# Velocity image
char	selimage[SZ_FNAME]	# Select image
char	britimage[SZ_FNAME]	# Brightness image
short	u[SWTH_HIGH], dat
int	date, hour, minute, seconds, i, j, num, lrs
pointer	velim, selim, britim, velsrp, selsrp, britsrp

int	read()
pointer	immap(), impl2s()
errchk	immap, impl2s

begin
	# Calculate the time.  Assemble the output image names.
	hour = int(VT_HTIME(hs)/3600)
	minute = int((VT_HTIME(hs) - hour * 3600)/60)
	seconds = int(VT_HTIME(hs) - hour * 3600 - minute * 60)
	if (brief) {
	    call sprintf (velimage[1], SZ_FNAME, "v%03d")
	        call pargi (filenumber)
	    call sprintf (selimage[1], SZ_FNAME, "s%03d")
	        call pargi (filenumber)
	    call sprintf (britimage[1], SZ_FNAME, "b%03d")
	        call pargi (filenumber)
	} else {
	    call sprintf (velimage[1], SZ_FNAME, "v%02d_%02d%02d_%03d")
	        call pargi (VT_HDAY(hs)) # day of month
	        call pargi (hour)
	        call pargi (minute)
	        call pargi (filenumber)
	    call sprintf (selimage[1], SZ_FNAME, "s%02d_%02d%02d_%03d")
	        call pargi (VT_HDAY(hs)) # day of month
	        call pargi (hour)
	        call pargi (minute)
	        call pargi (filenumber)
	    call sprintf (britimage[1], SZ_FNAME, "b%02d_%02d%02d_%03d")
	        call pargi (VT_HDAY(hs)) # day of month
	        call pargi (hour)
	        call pargi (minute)
	        call pargi (filenumber)
	}
	if (select) {
	    selim = immap (selimage, NEW_IMAGE, 0)
	    IM_NDIM(selim) = 2
	    IM_LEN(selim,1) = SWTH_HIGH/2
	    IM_LEN(selim,2) = VT_HNUMCOLS(hs)
	    IM_PIXTYPE(selim) = TY_SHORT
	    call imaddi (selim, "obs_time", VT_HTIME(hs))
	    date = VT_HMONTH(hs) * 10000 + VT_HDAY(hs) * 100 + VT_HYEAR(hs)
	    call imaddi (selim, "obs_date", date )
	    call imaddi (selim, "wv_lngth", VT_HWVLNGTH(hs))
	    call imaddi (selim, "obs_type", VT_HOBSTYPE(hs))
	    call imaddi (selim, "av_intns", VT_HAVINTENS(hs))
	    call imaddi (selim, "num_cols", VT_HNUMCOLS(hs))
	    call imaddi (selim, "intg/pix", VT_HINTGPIX(hs))
	    call imaddi (selim, "rep_time", VT_HREPTIME(hs))
	}
	if (bright) {
	    britim = immap (britimage, NEW_IMAGE, 0)
	    IM_NDIM(britim) = 2
	    IM_LEN(britim,1) = SWTH_HIGH/2
	    IM_LEN(britim,2) = VT_HNUMCOLS(hs)
	    IM_PIXTYPE(britim) = TY_SHORT
	    call imaddi (britim, "obs_time", VT_HTIME(hs))
	    date = VT_HMONTH(hs) * 10000 + VT_HDAY(hs) * 100 + VT_HYEAR(hs)
	    call imaddi (britim, "obs_date", date )
	    call imaddi (britim, "wv_lngth", VT_HWVLNGTH(hs))
	    call imaddi (britim, "obs_type", VT_HOBSTYPE(hs))
	    call imaddi (britim, "av_intns", VT_HAVINTENS(hs))
	    call imaddi (britim, "num_cols", VT_HNUMCOLS(hs))
	    call imaddi (britim, "intg/pix", VT_HINTGPIX(hs))
	    call imaddi (britim, "rep_time", VT_HREPTIME(hs))
	}
	if (velocity) {
	    velim = immap (velimage, NEW_IMAGE, 0)
	    IM_NDIM(velim) = 2
	    IM_LEN(velim,1) = SWTH_HIGH/2
	    IM_LEN(velim,2) = VT_HNUMCOLS(hs)
	    IM_PIXTYPE(velim) = TY_SHORT
	    call imaddi (velim, "obs_time", VT_HTIME(hs))
	    date = VT_HMONTH(hs) * 10000 + VT_HDAY(hs) * 100 + VT_HYEAR(hs)
	    call imaddi (velim, "obs_date", date )
	    call imaddi (velim, "wv_lngth", VT_HWVLNGTH(hs))
	    call imaddi (velim, "obs_type", VT_HOBSTYPE(hs))
	    call imaddi (velim, "av_intns", VT_HAVINTENS(hs))
	    call imaddi (velim, "num_cols", VT_HNUMCOLS(hs))
	    call imaddi (velim, "intg/pix", VT_HINTGPIX(hs))
	    call imaddi (velim, "rep_time", VT_HREPTIME(hs))
	}

	do j = 1, VT_HNUMCOLS(hs) {
	    if (select)
	        selsrp = impl2s (selim, j)
	    if (bright)
	        britsrp = impl2s (britim, j)
	    if (velocity)
	        velsrp = impl2s (velim, j)

	    iferr (num = read (inputfd, u, SWTH_HIGH*SZB_SHORT/SZB_CHAR)) {
		call fseti (inputfd, F_VALIDATE, lrs*SZB_SHORT/SZB_CHAR)
		call eprintf ("Error on tape read.\n")
		num = read (inputfd, u, SWTH_HIGH*SZB_SHORT/SZB_CHAR)
	    }
	    lrs = num
	    if (num < SWTH_HIGH*SZB_SHORT/SZB_CHAR)
	        call error (0, "eof encountered when reading file")
	    if (BYTE_SWAP2 == YES)
	        call bswap2 (u, 1, u, 1, SWTH_HIGH * SZB_SHORT)

	    # Unpack the data into the three arrays.
	    do i = 257, 512 {
		if (select) {
		    dat = u[i]/16
		    if (u[i] < 0)
			dat = dat - 1
		    Mems[selsrp+i-257] = dat
		}
		if (bright)
		    Mems[britsrp+i-257] = and(int(u[i]),17B)*16
	    }

	    do i = 1, 256 {
		if (velocity) {
		    dat = u[i]/16
		    if (u[i] < 0)
			dat = dat - 1
		    Mems[velsrp+i-1] = dat
		}
		if (bright)
		    Mems[britsrp+i-1] = Mems[britsrp+i-1]+and(int(u[i]),17B)
	    }
	}

	# Unmap images.
	if (select)
	    call imunmap (selim)
	if (velocity)
	    call imunmap (velim)
	if (bright)
	    call imunmap (britim)
end
