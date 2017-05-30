#define FUNCTIONS_H
char	*xmalloc(), *xcalloc();
double	**array ();
void	free_array ();
struct	token gettok();
struct	item item_copy(), pop();
void 	print_item(), print_stack(), print_top(),
	mark_stack(), reset_stack();
void    del_item(), store_name();
void 	dm_mul (), mm_mul (), tm_mul(), mt_mul (), tt_mul(), mm_plus(),
	mm_minus(), m_trn(), m_copy(), ident(), eigv(), m_shft(),
	m_hes(), m_qr(), m_qrp(), m_solve(), t_solve(), 
	m_svd(), m_bidiag(), chase (), sortsvs (), m_ginv(),
	control ();
double	m_rnorm(), m_cnorm(), m_norm();

/* Actual calculater functions */

void    alg1(), alg2(), alg3(), alg4(),
	plus(), minus(), times(), divide(), ev(), ginv(),
	print(), drop(), trn(), hes(), 
	get_mat(), dupl(), over(), swap(), svd(),
	quit(), qr(), qrp(), solve(), inverse(),
	idn(), fnrm(), store(), cont1(), cont2(),
        usin(), ucos(), utan(), uexp(), uatan(), uln(), ulog(),
        flops(), noflops(), usqrt(), uchs();
