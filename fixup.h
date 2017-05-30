#define fixup(r,s,mx)					\
	     ( (s != 0.0) ?				\
	       ( mx = max(fabs(r),fabs(s)),		\
		 r /= mx,				\
                 s /= mx,				\
	         mx = sqrt (r*r + s*s),			\
 	         r /= mx,		 	 	\
		 s /= mx ) :				\
               0                                        \
 	      )
