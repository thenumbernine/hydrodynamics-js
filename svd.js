/*******************************************************************************
 * Singular value decomposition program, svd, from "Numerical Recipes in C"
 * (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
 * and B.P. Flannery
 * *******************************************************************************/

function SIGN(a,b) { return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a)); }

function pythag(a, b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	var absa=Math.abs(a);
	var absb=Math.abs(b);
	if (absa > absb) return absa*Math.sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*Math.sqrt(1.0+(absa/absb)*(absa/absb)));
}

/******************************************************************************/
function svd(a, u, w, v)
/*******************************************************************************
 * Given a matrix a[0..m-1][0..n-1], this routine computes its singular value
 * decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
 * matrix of singular values W is output as a vector w[1..n].  The matrix V (not
 * the transpose VT) is output as v[1..n][1..n].
 * *******************************************************************************/
{
	var flag,i,its,j,jj,k,l,nm;
	var anorm,c,f,g,h,s,scale,x,y,z,rv1;

	var m = a.length;
	var n = a[0].length;
	
	u.length = m;
	for (var i = 0; i < m; ++i) {
		if (u[i] === undefined) u[i] = [];
		u[i].length = n;
		for (var j = 0; j < n; ++j) {
			u[i][j] = a[i][j];
		}
	}
	v.length = n;
	for (var i = 0; i < n; ++i) {
		if (v[i] === undefined) v[i] = [];
		v[i].length = n;
	}
	
	w.length = n; 

	rv1 = new Array(n);
	g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += Math.abs(u[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}
				f=u[i][i];
				g = -SIGN(Math.sqrt(s),f);
				h=f*g-s;
				u[i][i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
					f=s/h;
					for (k=i;k<m;k++) u[k][j] += f*u[k][i];
				}
				for (k=i;k<m;k++) u[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i < m && i != n-1) {
			for (k=l;k<n;k++) scale += Math.abs(u[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}
				f=u[i][l];
				g = -SIGN(Math.sqrt(s),f);
				h=f*g-s;
				u[i][l]=f-g;
				for (k=l;k<n;k++) rv1[k]=u[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += u[j][k]*u[i][k];
					for (k=l;k<n;k++) {
						u[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<n;k++) u[i][k] *= scale;
			}
		}
		var dmaxarg2 = Math.abs(w[i]) + Math.abs(rv1[i]);
		anorm = Math.max(anorm, dmaxarg2);
	}
	for (i=n-1;i>=0;i--) { /* Accumulation of right-hand transformations. */
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++) /* Double division to avoid possible underflow. */
					v[j][i]=(u[i][j] / u[i][l]) / g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=(m<n?m:n)-1;i>=0;i--) { /* Accumulation of left-hand transformations. */
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) u[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
				f=(s/u[i][i])*g;
				for (k=i;k<m;k++) u[k][j] += f*u[k][i];
			}
			for (j=i;j<m;j++) u[j][i] *= g;
		} else for (j=i;j<m;j++) u[j][i]=0.0;
		++u[i][i];
	}
	for (k=n-1;k>=0;k--) { /* Diagonalization of the bidiagonal form. */
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) { /* Test for splitting. */
				nm=l-1; /* Note that rv1[1] is always zero. */
				if ((Math.abs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((Math.abs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; /* Cancellation of rv1[l], if l > 1. */
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((Math.abs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u[j][nm];
						z=u[j][i];
						u[j][nm]=y*c+z*s;
						u[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) throw "no convergence in 30 svd iterations";
			x=w[l]; /* Shift from bottom 2-by-2 minor. */
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; /* Next QR transformation: */
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

