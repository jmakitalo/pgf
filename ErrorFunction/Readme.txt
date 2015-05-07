
				     MATLAB toolbox
		________________________________________________________

			   Error function of complex numbers
		________________________________________________________


** Contents

	1. Overview
	2. Requirements
	3. Installation
	4. Copyright
	5. Warranty
	6. History
	7. Download
	8. Trademarks

** Publisher

	Marcel Leutenegger		marcel.leutenegger@a3.epfl.ch
	EPFL STI SMT LOB
	BM 4.245			Phone:	+41 21 693 78 21
	Station 17
	CH-1015 Lausanne



1. Overview

	This package provides improved implementations of the error function for
	MATLAB. It ships a MEX-file for calculating the	error function of real-
	valued numbers 5-6x faster than with the default MATLAB implementation.
	A second MEX-file and/or a companion M-file enhances the calculation of
	the error function for complex numbers. See "erfz.pdf" for implementa-
	tion details.


2. Requirements

	• An Intel Pentium II compatible computer or newer.
	• MATLAB 6.0 or newer running.


3. Installation

	Unpack the archive in a folder that is part of the MATLAB path. The error
	functions should reside in a '@double' folder to avoid potential data type
	conflicts.


4. Copyright

	This software is published as freeware. The author reserves the right to
	modify any of the contained files.

	You are allowed to distribute the functions as long as you deliver for free the
	entire package.

		Path		Files

		/		Readme.txt
		@double/	erf.dll		Replaces MATLAB's error function
				erfz.dll	Error function of complex numbers
				erfz.m		Help and companion function
				erfz.pdf	Implementation details


5. Warranty

	Any warranty is strictly refused and you cannot anticipate any financial or
	technical support in case of malfunction or damage.

	Feedback and comments are welcome. I will try to track reported problems and
	fix bugs.


6. History

   • January 14, 2008
	Initial relase.


7. Download

	Optimized MATLAB functions are available online at (subject to change):

		   http://ioalinux1.epfl.ch/~mleutene/MATLABToolbox/


	Summaries are also published at MATLAB central:

			http://www.mathworks.com/matlabcentral/


8. Trademarks

	MATLAB is a registered trademark of The MathWorks, Inc. Pentium is a
	registered trademark of Intel Corporation. Other product or brand names
	are trademarks or registered trademarks of their respective holders.

		________________________________________________________

			    Site map • EPFL © 2008, Lausanne
			       Webmaster • 14 Janaury 2008
