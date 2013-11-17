/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured
*/

var panel;
var canvas;
var waveVtxBuf, waveStateBuf;
var xmin = 0;
var xmax = 100; 
var ymin = 0;
var ymax = 100;
var gridstep = 10;
var mouse;

//interface directions
var dirs = [[1,0], [0,1]];

function isnan(x) { x = Number(x); return x != x; }

function mat33invert(out, a) {
	var det = a[0][0] * a[1][1] * a[2][2]
			+ a[1][0] * a[2][1] * a[0][2]
			+ a[2][0] * a[0][1] * a[1][2]
			- a[2][0] * a[1][1] * a[0][2]
			- a[1][0] * a[0][1] * a[2][2]
			- a[0][0] * a[2][1] * a[1][2];
	if (det == 0) {
		for (var j = 0; j < 3; ++j) {
			for (var i = 0; i < 3; ++i) {
				console.log('a('+i+','+j+') = '+a[j][i]);
			}
		}
		throw 'singular!';
	}
	var invDet = 1 / det;
	for (var j = 0; j < 3; ++j) {
		var j1 = (j + 1) % 3;
		var j2 = (j + 2) % 3;
		for (var i = 0; i < 3; ++i) {
			var i1 = (i + 1) % 3;
			var i2 = (i + 2) % 3;
			out[i][j] = invDet * (a[j1][i1] * a[j2][i2] - a[j1][i2] * a[j2][i1]);
		}
	}
}

/*
output:
matrix
eigenvalues
eigenvectors

input:
velocity
hTotal
speedOfSound
*/
function buildEigenstate(offset, matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocityX, velocityY, hTotal, gamma, normalX, normalY) {

	if ((hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)) < 0) throw 'sqrt error';

	//calculate matrix & eigenvalues & vectors at interface from state at interface
	var speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)));
	var tangentX = -normalY;
	var tangentY = normalX;
	var velocityN = velocityX * normalX + velocityY * normalY;
	var velocityT = velocityX * tangentX + velocityY * tangentY;
	var velocitySq = velocityX * velocityX + velocityY * velocityY;	
	
	//eigenvalues: min, mid, max
	eigenvalues[0 + 4 * offset] = velocityN - speedOfSound;
	eigenvalues[1 + 4 * offset] = velocityN;
	eigenvalues[2 + 4 * offset] = velocityN;
	eigenvalues[3 + 4 * offset] = velocityN + speedOfSound;

	//mid/tangent eigenvector, velocity components:
	//http://www.cfdbooks.com/cfdcodes has a speed of sound scale ...
	//"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro does not
	//http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_7.pdf also does not
	//http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf also does not

	//mid/tangent eigenvector, energy component:
	//http://www.cfdbooks.com/cfdcodes has tangent velocity times speed of sound (he sure likes scaling things by the speed of sound)
	//"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro has the tangent velocity
	//http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_7.pdf has the tangent velocity squared
	//http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf has tangent velocity 

	//I'm going with http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf
	//...but it looks like their tangent vector is <y,-x>, (-90' rotation)
	// whereas I'm going to use <-y,x> (+90' rotation)

	//min eigenvector
	eigenvectors[0 + 4 * (0 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (0 + 4 * offset)] = velocityX - speedOfSound * normalX;
	eigenvectors[2 + 4 * (0 + 4 * offset)] = velocityY - speedOfSound * normalY;
	eigenvectors[3 + 4 * (0 + 4 * offset)] = hTotal - speedOfSound * velocityN;
	//mid eigenvector (normal)
	eigenvectors[0 + 4 * (1 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (1 + 4 * offset)] = velocityX;
	eigenvectors[2 + 4 * (1 + 4 * offset)] = velocityY;
	eigenvectors[3 + 4 * (1 + 4 * offset)] = .5 * velocitySq;
	//mid eigenvector (tangent)
	eigenvectors[0 + 4 * (2 + 4 * offset)] = 0;
	eigenvectors[1 + 4 * (2 + 4 * offset)] = tangentX;
	eigenvectors[2 + 4 * (2 + 4 * offset)] = tangentY;
	eigenvectors[3 + 4 * (2 + 4 * offset)] = velocityT;
	//max eigenvector
	eigenvectors[0 + 4 * (3 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (3 + 4 * offset)] = velocityX + speedOfSound * normalX;
	eigenvectors[2 + 4 * (3 + 4 * offset)] = velocityY + speedOfSound * normalY;
	eigenvectors[3 + 4 * (3 + 4 * offset)] = hTotal + speedOfSound * velocityN;
	
	//calculate eigenvector inverses ... 
	//min row
	eigenvectorsInverse[0 + 4 * (0 + 4 * offset)] = (.5 * (gamma - 1) * velocitySq + speedOfSound * velocityN) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (1 + 4 * offset)] = -(normalX * speedOfSound + (gamma - 1) * velocityX) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (2 + 4 * offset)] = -(normalY * speedOfSound + (gamma - 1) * velocityY) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (3 + 4 * offset)] = (gamma - 1) / (2 * speedOfSound * speedOfSound);
	//mid normal row
	eigenvectorsInverse[1 + 4 * (0 + 4 * offset)] = 1 - .5 * (gamma - 1) * velocitySq / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (1 + 4 * offset)] = (gamma - 1) * velocityX / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (2 + 4 * offset)] = (gamma - 1) * velocityY / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (3 + 4 * offset)] = -(gamma - 1) / (speedOfSound * speedOfSound);
	//mid tangent row
	eigenvectorsInverse[2 + 4 * (0 + 4 * offset)] = -velocityT; 
	eigenvectorsInverse[2 + 4 * (1 + 4 * offset)] = tangentX;
	eigenvectorsInverse[2 + 4 * (2 + 4 * offset)] = tangentY;
	eigenvectorsInverse[2 + 4 * (3 + 4 * offset)] = 0;
	//max row
	eigenvectorsInverse[3 + 4 * (0 + 4 * offset)] = (.5 * (gamma - 1) * velocitySq - speedOfSound * velocityN) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (1 + 4 * offset)] = (normalX * speedOfSound - (gamma - 1) * velocityX) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (2 + 4 * offset)] = (normalY * speedOfSound - (gamma - 1) * velocityY) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (3 + 4 * offset)] = (gamma - 1) / (2 * speedOfSound * speedOfSound);

	//calculate matrix
	for (var i = 0; i < 4; ++i) {
		for (var j = 0; j < 4; ++j) {
			var s = 0;
			var identCheck = 0;
			for (var k = 0; k < 4; ++k) {
				s += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)] * eigenvalues[k + 4 * offset];
				identCheck += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)];
			}
			matrix[i + 4 * (j + 4 * offset)] = s;
			/** /
			var epsilon = 1e-5;
			if (Math.abs(identCheck - (i == j ? 1 : 0)) > epsilon) throw 'bad eigen basis';
			/**/
		}
	}
	/** /	
	function f32subset(a, o, s) {
		var d = new Float32Array(s);
		for (var i = 0; i < s; ++i) {
			d[i] = a[i+o];
		}
		return d;
	}
	console.log('offset',offset);
	console.log('velocity',velocityX,velocityY);
	console.log('hTotal',hTotal);
	console.log('gamma',gamma);
	console.log('normal',normalX,normalY);

	console.log('eigenvalues:',f32subset(eigenvalues, 4*offset, 4));
	console.log('eigenvectors:',f32subset(eigenvectors, 16*offset, 16));
	console.log('eigenvectors^-1:',f32subset(eigenvectorsInverse, 16*offset, 16));
	console.log('matrix:',f32subset(matrix, 16*offset, 16));
	console.log('e^-1 * e:',identCheck);
	throw 'here';
	/**/
}

var fluxMethods = {
	donorCell : function(r) { return 0; },
	laxWendroff : function(r) { return 1; },
	
	beamWarming : function(r) { return r; },
	fromm : function(r) { return .5 * (1 + r); },

	//Wikipedia
	CHARM : function(r) { return Math.max(0, r*(3*r+1)/((r+1)*(r+1)) ); },
	HCUS : function(r) { return Math.max(0, 1.5 * (r + Math.abs(r)) / (r + 2) ); },
	HQUICK : function(r) { return Math.max(0, 2 * (r + Math.abs(r)) / (r + 3) ); },
	Koren : function(r) { return Math.max(0, Math.min(2*r, (1 + 2*r)/3 ,2) ); },
	minmod : function(r) { return Math.max(0, Math.min(r,1) ); },
	Oshker : function(r) { return Math.max(0, Math.min(r,1.5) ); },	//replace 1.5 with 1 <= beta <= 2	
	ospre : function(r) { return .5 * (r*r + r) / (r*r + r + 1); },
	smart : function(r) { return Math.max(0, Math.min(2 * r, .25 + .75 * r, 4)); },
	Sweby : function(r) { return Math.max(0, Math.min(1.5 * r, 1), Math.min(r, 1.5)); },	//replace 1.5 with 1 <= beta <= 2
	UMIST : function(r) { return Math.max(0, Math.min(2*r, .75 + .25*r, .25 + .75*r, 2)); },	
	vanAlbada1 : function(r) { return (r * r + r) / (r * r + 1); },
	vanAlbada2 : function(r) { return 2 * r / (r * r + 1); },
	
	vanLeer : function(r) { return (r + Math.abs(r)) / (1 + Math.abs(r)); },
	MC : function(r) { return Math.max(0, Math.min(2, .5 * (1 + r), 2 * r)); },
	superbee : function(r) { return Math.max(0,Math.min(1,2*r),Math.min(2,r)); }
};

var boundaryMethods = {
	periodic : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var state = 0; state < 4; ++state) {
				//top
				q[state + 4 * (i + nx * 0)] = q[state + 4 * (i + nx * (nx-4))];
				q[state + 4 * (i + nx * 1)] = q[state + 4 * (i + nx * (nx-3))];
				q[state + 4 * (i + nx * (nx-2))] = q[state + 4 * (i + nx * 2)];
				q[state + 4 * (i + nx * (nx-1))] = q[state + 4 * (i + nx * 3)];
				//left	
				q[state + 4 * (0 + nx * i)] = q[state + 4 * (nx-4 + nx * i)];
				q[state + 4 * (1 + nx * i)] = q[state + 4 * (nx-3 + nx * i)];
				q[state + 4 * (nx-2 + nx * i)] = q[state + 4 * (2 + nx * i)];
				q[state + 4 * (nx-1 + nx * i)] = q[state + 4 * (3 + nx * i)];
			}
		}
	},
	mirror : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			//top
			q[0 + 4 * (i + nx * (0))] = q[0 + 4 * (i + nx * (3))];
			q[0 + 4 * (i + nx * (1))] = q[0 + 4 * (i + nx * (2))];
			q[0 + 4 * (i + nx * (nx-2))] = q[0 + 4 * (i + nx * (nx-3))];
			q[0 + 4 * (i + nx * (nx-1))] = q[0 + 4 * (i + nx * (nx-4))];
			q[1 + 4 * (i + nx * (0))] = -q[1 + 4 * (i + nx * (3))];
			q[1 + 4 * (i + nx * (1))] = -q[1 + 4 * (i + nx * (2))];
			q[1 + 4 * (i + nx * (nx-2))] = -q[1 + 4 * (i + nx * (nx-3))];
			q[1 + 4 * (i + nx * (nx-1))] = -q[1 + 4 * (i + nx * (nx-4))];
			q[2 + 4 * (i + nx * (0))] = -q[2 + 4 * (i + nx * (3))];
			q[2 + 4 * (i + nx * (1))] = -q[2 + 4 * (i + nx * (2))];
			q[2 + 4 * (i + nx * (nx-2))] = -q[2 + 4 * (i + nx * (nx-3))];
			q[2 + 4 * (i + nx * (nx-1))] = -q[2 + 4 * (i + nx * (nx-4))];
			q[3 + 4 * (i + nx * (0))] = q[3 + 4 * (i + nx * (3))];
			q[3 + 4 * (i + nx * (1))] = q[3 + 4 * (i + nx * (2))];
			q[3 + 4 * (i + nx * (nx-2))] = q[3 + 4 * (i + nx * (nx-3))];
			q[3 + 4 * (i + nx * (nx-1))] = q[3 + 4 * (i + nx * (nx-4))];
			//left
			q[0 + 4 * (0 + nx * i)] = q[0 + 4 * (3 + nx * i)];
			q[0 + 4 * (1 + nx * i)] = q[0 + 4 * (2 + nx * i)];
			q[0 + 4 * (nx-2 + nx * i)] = q[0 + 4 * (nx-3 + nx * i)];
			q[0 + 4 * (nx-1 + nx * i)] = q[0 + 4 * (nx-4 + nx * i)];
			q[1 + 4 * (0 + nx * i)] = -q[1 + 4 * (3 + nx * i)];
			q[1 + 4 * (1 + nx * i)] = -q[1 + 4 * (2 + nx * i)];
			q[1 + 4 * (nx-2 + nx * i)] = -q[1 + 4 * (nx-3 + nx * i)];
			q[1 + 4 * (nx-1 + nx * i)] = -q[1 + 4 * (nx-4 + nx * i)];
			q[2 + 4 * (0 + nx * i)] = -q[2 + 4 * (3 + nx * i)];
			q[2 + 4 * (1 + nx * i)] = -q[2 + 4 * (2 + nx * i)];
			q[2 + 4 * (nx-2 + nx * i)] = -q[2 + 4 * (nx-3 + nx * i)];
			q[2 + 4 * (nx-1 + nx * i)] = -q[2 + 4 * (nx-4 + nx * i)];
			q[3 + 4 * (0 + nx * i)] = q[3 + 4 * (3 + nx * i)];
			q[3 + 4 * (1 + nx * i)] = q[3 + 4 * (2 + nx * i)];
			q[3 + 4 * (nx-2 + nx * i)] = q[3 + 4 * (nx-3 + nx * i)];
			q[3 + 4 * (nx-1 + nx * i)] = q[3 + 4 * (nx-4 + nx * i)];
		}
	},
	dirichlet : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var state = 0; state < 4; ++state) {
				q[state + 4 * (i + nx * (0))] = 0;
				q[state + 4 * (1 + nx * (1))] = 0;
				q[state + 4 * (i + nx * (nx-2))] = 0;
				q[state + 4 * (i + nx * (nx-1))] = 0;
				q[state + 4 * (0 + nx * i)] = 0;
				q[state + 4 * (1 + nx * (1))] = 0;
				q[state + 4 * (nx-2 + nx * i)] = 0;
				q[state + 4 * (nx-1 + nx * i)] = 0;
			}
		}
	},
	constant : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var state = 0; state < 4; ++state) {
				q[state + 4 * (i + nx * (0))] = q[state + 4 * (i + nx * (1))] = q[state + 4 * (i + nx * (2))];
				q[state + 4 * (i + nx * (nx-1))] = q[state + 4 * (i + nx * (nx-2))] = q[state + 4 * (i + nx * (nx-3))];
				q[state + 4 * (0 + nx * i)] = q[state + 4 * (1 + nx * i)] = q[state + 4 * (2 + nx * i)];
				q[state + 4 * (nx-1 + nx * i)] = q[state + 4 * (nx-2 + nx * i)] = q[state + 4 * (nx-3 + nx * i)];
			}
		}
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {
			var mindum = undefined;
			var qIndex = 0;
			for (var j = 0; j < this.nx; ++j) {
				for (var i = 0; i < this.nx; ++i) {
					var rho = this.q[0 + qIndex];
					var u = this.q[1 + qIndex] / rho;
					var v = this.q[2 + qIndex] / rho; 
					var energyTotal = this.q[3 + qIndex] / rho; 
					var energyKinematic = .5 * (u * u + v * v);
					var energyThermal = energyTotal - energyKinematic;
					var speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
					var dx = this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] - this.xi[0 + 2 * (i + (this.nx+1) * j)];
					var dy = this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] - this.xi[1 + 2 * (i + (this.nx+1) * j)];
					//http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
					//var dum = 1 / (speedOfSound + Math.abs(u) / dx + Math.abs(v) / dy); 
					var dum = dx / (speedOfSound + Math.abs(u));
					if (mindum === undefined || dum < mindum) mindum = dum;
					var dum = dy / (speedOfSound + Math.abs(v));
					if (mindum === undefined || dum < mindum) mindum = dum;
					qIndex += 4;
				}
			}
			return this.cfl * mindum;
		},
		advect : function(dt) {
			//assert(this.x.length == this.nx);
			//assert(this.xi.length == this.nx + 1);
			//assert(this.q.length == this.nx);
			//assert(this.ui.length == this.nx + 1);
	
			//get velocity at interfaces from state
			for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
				for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
					var iu = 2 * (i + (this.nx+1) * j);
					for (var side = 0; side < 2; ++side) {
						var iq1 = 4 * (i + this.nx * j);
						var iq2 = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
						this.ui[side + iu] = .5 * (
							this.q[1+side + iq1] / this.q[0 + iq1] + 
							this.q[1+side + iq2] / this.q[0 + iq2]);
					}
				}
			}
			//boundary zero
			for (var j = 0; j < this.nghost; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					for (var side = 0; side < 2; ++side) {
						//left boundary, left and top sides, zero vector
						this.ui[side + 2 * (j + (this.nx+1) * i)] = 0;
						//right boundary, left and top sides, zero vector
						this.ui[side + 2 * (this.nx-j + (this.nx+1) * i)] = 0;
						//top boundary, left and top sides, zero vector
						this.ui[side + 2 * (i + (this.nx+1) * j)] = 0;
						//bottom boundary, left and top sides, zero vector
						this.ui[side + 2 * (i + (this.nx+1) * (this.nx-j))] = 0;
					}
				}
			}

			//compute flux and advect for each state vector
			for (var state = 0; state < 4; ++state) {	//state
				//r_{i-1/2},{j-1/2} flux limiter
				for (var j = this.nghost; j < this.nx+this.nghost-3; ++j) {
					for (var i = this.nghost; i < this.nx+this.nghost-3; ++i) {
						for (var side = 0; side < 2; ++side) {
							//dq = q_i,j - q_{{i,j}-dirs[side]}
							var dq = this.q[state + 4 * (i + this.nx * j)]
								- this.q[state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))];
							if (Math.abs(dq) > 0) {
								if (this.ui[side + 2 * (i + (this.nx+1) * j)] >= 0) {
									this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))]
											- this.q[state + 4 * (i - 2*dirs[side][0] + this.nx * (j - 2*dirs[side][1]))]) / dq;
								} else {
									this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[state + 4 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]))]
											- this.q[state + 4 * (i + this.nx * j)]) / dq;
									
								}
							} else {
								this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							}
						}
					}
				}
			
				//now for ghost boundaries
				for (var j = 0; j < this.nghost; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.r[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.r[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.r[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						}
					}
				}


				var sideLengths = [];	//[2]
				//construct flux:
				for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						var dx = this.x[0 + 2 * (i + this.nx * j)] - this.x[0 + 2 * (i-1 + this.nx * j)];
						var dy = this.x[1 + 2 * (i + this.nx * j)] - this.x[1 + 2 * (i + this.nx * (j-1))];
						var volume = dx * dy;
						sideLengths[0] = dx;
						sideLengths[1] = dy;
					
						//flux calculation 
						for (var side = 0; side < 2; ++side) {
							var ir = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							//apply limiter
							var phi = this.fluxMethod(this.r[ir]);
							var iu = side + 2 * (i + (this.nx+1) * j);
							var iflux = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							var iqL = state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
							var iqR = state + 4 * (i + this.nx * j);
							if (this.ui[iu] >= 0) {
								this.flux[iflux] = this.ui[iu] * this.q[iqL];
							} else {
								this.flux[iflux] = this.ui[iu] * this.q[iqR];
							}
							var delta = phi * (this.q[iqR] - this.q[iqL]);
							this.flux[iflux] += delta * .5
								* Math.abs(this.ui[iu])
								* (1 - Math.abs(this.ui[iu] * dt * volume / (sideLengths[side] * sideLengths[side])));
						}
					}
				}
				//now for ghost boundaries
				for (var j = 0; j < this.nghost-1; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.flux[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.flux[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						}
					}
				}

				//update cells
				for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
						for (var side = 0; side < 2; ++side) {
							var ixR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
							var ixL = side + 2 * (i + (this.nx+1) * j);
							var dx = this.xi[ixR] - this.xi[ixL];
							var ifluxR = state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
							var ifluxL = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							var df = this.flux[ifluxR] - this.flux[ifluxL];
							this.q[state + 4 * (i + this.nx * j)] -= dt * df / dx;
						}
					}
				}
			}
		}
	},
	Riemann : {
		initStep : function() {
			var mindum = undefined;
			for (var iy = 1; iy < this.nx; ++iy) {
				for (var ix = 1; ix < this.nx; ++ix) {
					for (var side = 0; side < 2; ++side) {
						var qIndexL = 4 * (ix - dirs[side][0] + this.nx * (iy - dirs[side][1]));
						var densityL = this.q[0 + qIndexL];
						var velocityXL = this.q[1 + qIndexL] / densityL;
						var velocityYL = this.q[2 + qIndexL] / densityL;
						var energyTotalL = this.q[3 + qIndexL] / densityL;
						var energyKinematicL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
						var energyThermalL = energyTotalL - energyKinematicL;
						var pressureL = (this.gamma - 1) * densityL * energyThermalL;
						var speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
						var hTotalL = energyTotalL + pressureL / densityL;
						var roeWeightL = Math.sqrt(densityL);
						
						var qIndexR = 4 * (ix + this.nx * iy);
						var densityR = this.q[0 + qIndexR];
						var velocityXR = this.q[1 + qIndexR] / densityR;
						var velocityYR = this.q[2 + qIndexR] / densityR;
						var energyTotalR = this.q[3 + qIndexR] / densityR;
						var energyKinematicR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
						var energyThermalR = energyTotalR - energyKinematicR;
						var pressureR = (this.gamma - 1) * densityR * energyThermalR;
						var speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
						var hTotalR = energyTotalR + pressureR / densityR;
						var roeWeightR = Math.sqrt(densityR);

						var denom = roeWeightL + roeWeightR;
						var velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
						var velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
						var hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
						if (isnan(velocityX)) throw 'nan';
						if (isnan(velocityY)) throw 'nan';
						if (isnan(hTotal)) throw 'nan';

						if ((hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)) < 0) throw 'cannot calculate speed of sound';
						
						buildEigenstate(
							 //index into interface element.  
							 //from there you'll have to scale by cell size.  
							 //Thus manually recreating the automatic storage of C structures. 
							 //JavaScript, why can't you be more like LuaJIT? 
							 side + 2 * (ix + (this.nx+1) * iy),	
							 this.interfaceMatrix,	//dim^2 = 16
							 this.interfaceEigenvalues,	//dim = 4
							 this.interfaceEigenvectors,	//dim^2 = 16
							 this.interfaceEigenvectorsInverse,	//dim^2 = 16
							 velocityX, velocityY, hTotal, this.gamma,
							 dirs[side][0], dirs[side][1]);
						
						if (isnan(this.interfaceEigenvalues[0 + 4 * (side + 2 * (ix + (this.nx+1) * iy))])) throw 'nan';
						if (isnan(this.interfaceEigenvalues[1 + 4 * (side + 2 * (ix + (this.nx+1) * iy))])) throw 'nan';
						if (isnan(this.interfaceEigenvalues[2 + 4 * (side + 2 * (ix + (this.nx+1) * iy))])) throw 'nan';
						if (isnan(this.interfaceEigenvalues[3 + 4 * (side + 2 * (ix + (this.nx+1) * iy))])) throw 'nan';

						var maxLambda = Math.max(0, 
							this.interfaceEigenvalues[0+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[1+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[2+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[3+4*(side+2*(ix+(this.nx+1)*iy))]);
						var minLambda = Math.min(0, 
							this.interfaceEigenvalues[0+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[1+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[2+4*(side+2*(ix+(this.nx+1)*iy))],
							this.interfaceEigenvalues[3+4*(side+2*(ix+(this.nx+1)*iy))]);
						var dx = this.xi[side + 2 * (ix+dirs[side][0] + (this.nx+1) * (iy+dirs[side][1]))] 
							- this.xi[side + 2 * (ix + (this.nx+1) * iy)];
						var dum = dx / (maxLambda - minLambda);
						if (mindum === undefined || dum < mindum) mindum = dum;
					}
				}
			}
			if (isnan(mindum)) throw 'nan';
			return this.cfl * mindum;
	
		},
		advect : function(dt) {
			for (var iy = 1; iy < this.nx; ++iy) {
				for (var ix = 1; ix < this.nx; ++ix) {
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							//find state change across interface in the basis of the eigenspace at the interface
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
									//reproject into interface eigenspace
								sum += this.interfaceEigenvectorsInverse[state + 4 * (k + 4 * (side + 2 * (ix + (this.nx+1) * iy)))]
									//flux difference
									* (this.q[k + 4 * (ix + this.nx * iy)] 
										- this.q[k + 4 * (ix - dirs[side][0] + this.nx * (iy - dirs[side][1]))])
							}
							this.interfaceDeltaQTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] = sum;
						}
					}
				}
			}

			for (var iy = this.nghost; iy < this.nx + this.nghost - 3; ++iy) {
				for (var ix = this.nghost; ix < this.nx + this.nghost - 3; ++ix) {
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							var interfaceDeltaQTilde = this.interfaceDeltaQTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))];
							if (Math.abs(interfaceDeltaQTilde) > 0) {
								if (this.interfaceEigenvalues[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] > 0) {
									this.rTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] = 
										this.interfaceDeltaQTilde[state + 4 * (side + 2 * (ix - dirs[side][0] + (this.nx+1) * (iy - dirs[side][1])))]
										/ interfaceDeltaQTilde;
								} else {
									this.rTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] = 
										this.interfaceDeltaQTilde[state + 4 * (side + 2 * (ix + dirs[side][0] + (this.nx+1) * (iy + dirs[side][1])))]
										/ interfaceDeltaQTilde;
								}
							}
						}
					}
				}
			}

			//..and keep the boundary rTilde's zero	
			for (var j = 0; j < this.nghost; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.rTilde[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.rTilde[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						}
					}
				}
			}
		
			var fluxAvg = [];	//4
			var fluxTilde = [];	//4
			var sideLengths = [];	//2
			
			//transform cell q's into cell qTilde's (eigenspace)
			for (var iy = 1; iy < this.nx; ++iy) {
				for (var ix = 1; ix < this.nx; ++ix) {
					var dx = this.x[0 + 2 * (i + this.nx * j)] - this.x[0 + 2 * (i-1 + this.nx * j)];
					var dy = this.x[1 + 2 * (i + this.nx * j)] - this.x[1 + 2 * (i + this.nx * (j-1))];
					var volume = dx * dy;
					sideLengths[0] = dx;
					sideLengths[1] = dy;
				
					for (var side = 0; side < 2; ++side) {
						
						//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
						//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
						for (var state = 0; state < 4; ++state) {
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
								sum += this.interfaceMatrix[state + 4 * (k + 4 * (side + 2 * (ix + (this.nx+1) * iy)))]
									* (this.q[k + 4 * (ix - dirs[side][0] + this.nx * (iy - dirs[side][1]))]
										+ this.q[k + 4 * (ix + this.nx * iy)]);
							}
							fluxAvg[state] = .5 * sum;
						}

						//calculate flux
						for (var state = 0; state < 4; ++state) {
							var theta = 0;
							if (this.interfaceEigenvalues[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] >= 0) {
								theta = 1;
							} else {
								theta = -1;
							}
						
							var phi = this.fluxMethod(this.rTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))]);
							if (isnan(phi)) throw 'nan';

							var epsilon = this.interfaceEigenvalues[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] * dt * volume / (sideLengths[side] * sideLengths[side]);
							if (isnan(epsilon)) throw 'nan';

							var deltaFluxTilde = this.interfaceEigenvalues[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))]
								* this.interfaceDeltaQTilde[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))];
							if (isnan(deltaFluxTilde)) throw 'nan';

							fluxTilde[state] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
							if (isnan(fluxTilde[state])) throw 'nan';
						}
					
						//reproject fluxTilde back into q
						for (var state = 0; state < 4; ++state) {
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
								sum += fluxTilde[k]
									* this.interfaceEigenvectors[state + 4 * (k + 4 * (side + 2 * (ix + (this.nx+1) * iy)))];
							}
							if (isnan(sum)) throw 'nan';
							this.flux[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))] = fluxAvg[state] + sum;
							if (isnan(this.flux[state + 4 * (side + 2 * (ix + (this.nx+1) * iy))])) throw 'nan';
						}
					}
				}
			}
		
			//zero boundary flux
			//..and keep the boundary r's zero	
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * 0))] = 0;
						this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * this.nx))] = 0;
						this.flux[state + 4 * (side + 2 * (0 + (this.nx+1) * i))] = 0;
						this.flux[state + 4 * (side + 2 * (this.nx + (this.nx+1) * i))] = 0;
					}
				}
			}

			//update cells
			for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
					for (var side = 0; side < 2; ++side) {
						var xiIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
						var xiIndex = side + 2 * (i + (this.nx+1) * j);
						var dx = this.xi[xiIndexR] - this.xi[xiIndex];
						for (var state = 0; state < 4; ++state) {
							this.q[state + 4 * (i + this.nx * j)] -= 
								(this.flux[state + 4 * (side + 2 * (i+dirs[side][0] + (this.nx+1) * (j+dirs[side][0])))]
									- this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * j))]) * dt / dx;
							if (isnan(this.q[state + 4 * (i + this.nx * j)])) throw 'nan';	
						}
					}
				}
			}
		}
	}
};

var HydroState = makeClass({ 
	init : function() {
		this.nx = 200;
		this.cfl =.5;
		this.gamma = 7/5;
		
		//x_i,j,dim: cell positions
		//0 <= i < this.nx
		this.x = new Float32Array(this.nx * this.nx * 2);
		var e = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.x[e] = xmin + (xmax - xmin) * i / (this.nx-1); ++e;
				this.x[e] = ymin + (ymax - ymin) * j / (this.nx-1); ++e;
			}
		}
		
		//x_{i-1/2},{j-1/2},dim: interface positions
		//0 <= i,j < this.nx+1
		//0 <= dim < 2
		this.xi = new Float32Array((this.nx+1) * (this.nx+1) * 2);
		var e = 0;
		for (var j = 0; j <= this.nx; ++j) {
			var J = Math.min(j, this.nx-1);
			for (var i = 0; i <= this.nx; ++i) {
				var I = Math.min(i, this.nx-1);
				if (i == 0) {
					this.xi[e] = 2 * this.x[0 + 2 * (1 + this.nx * J)] - this.x[0 + 2 * (2 + this.nx * J)];
				} else if (i == this.nx) {
					this.xi[e] = 2 * this.x[0 + 2 * (this.nx-1 + this.nx * J)] - this.x[0 + 2 * (this.nx-2 + this.nx * J)];
				} else {
					this.xi[e] = .5*(this.x[0 + 2 * (i + this.nx * J)] + this.x[0 + 2 * (i-1 + this.nx * J)]);
				}
				++e;
				if (j == 0) {
					this.xi[e] = 2 * this.x[1 + 2 * (I + this.nx * 1)] - this.x[1 + 2 * (I + this.nx * 2)];
				} else if (j == this.nx) {
					this.xi[e] = 2 * this.x[1 + 2 * (I + this.nx * (this.nx-1))] - this.x[1 + 2 * (I + this.nx * (this.nx-2))];
				} else {
					this.xi[e] = .5*(this.x[1 + 2 * (I + this.nx * j)] + this.x[1 + 2 * (I + this.nx * (j-1))]);
				}
				++e;
			}
		}

		//q_i,j,state: state vector, stored as q[state + 4 * (j + this.nx * i)]
		//q_i,j,0: density: rho
		//q_i,j,1: momentum: rho * u
		//q_i,j,2: momentum: rho * v
		//q_i,j,3: work: rho * e
		this.q = new Float32Array(this.nx * this.nx * 4);
		var e = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				for (var state = 0; state < 4; ++state) {
					this.q[e] = 0; ++e;
				}
			}
		}

		this.resetSod();
		
		//p_i,j: pressure
		this.pressure = new Float32Array(this.nx * this.nx);
		var e = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.pressure[e] = 0; ++e;
			}
		}

		//TODO it is tempting to merge r, f, and ui into an edge structure
		//and associate them with the nodes on either side of them,
		//but then I would lose out on the 2nd-order contributions to the flux limiter.


		//f_{i-1/2},{j-1/2},side,state: cell flux
		this.flux = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4);
		var e = 0;
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.flux[e] = 0; ++e;
					}
				}
			}
		}
	

		//used for Burgers
		
		
		//r_{i-1/2},{j-1/2},side,state	
		this.r = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4);;
		var e = 0;
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.r[e] = 0; ++e;
					}
				}
			}
		}
		
	
		//only used with Burger's eqn advection code
		//u_{i-1/2},{j-1/2},dim: interface velocity
		this.ui = new Float32Array((this.nx+1) * (this.nx+1) * 2);
		var e = 0;
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					this.ui[e] = 0; ++e;
				}
			}
		}


		//used for Riemann
	

		//a_{i-1/2},{j-1/2},side,state,state
		this.interfaceMatrix = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		this.interfaceEigenvalues = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4);
		this.interfaceEigenvectors = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		this.interfaceEigenvectorsInverse = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var stateJ = 0; stateJ < 4; ++stateJ) {
						for (var stateI = 0; stateI < 4; ++stateI) {
							//initialize to identity matrix
							this.interfaceMatrix[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
							this.interfaceEigenvectors[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
							this.interfaceEigenvectorsInverse[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
						}
						this.interfaceEigenvalues[stateJ + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

		//qiTilde_{i-1/2},{j-1/2},side,state	
		this.interfaceDeltaQTilde = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4);
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

		//rTilde_{i-1/2},{j-1/2},side,state
		this.rTilde = new Float32Array((this.nx+1) * (this.nx+1) * 2 * 4);
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethod = boundaryMethods.mirror;
		this.fluxMethod = fluxMethods.superbee;
		this.advectMethod = advectMethods.Burgers;
	},
	resetSod : function() {
		var xIndex = 0;
		var qIndex = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.q[0 + qIndex] = (this.x[0 + xIndex] < 30 && this.x[1 + xIndex] < 30) ? 1 : .1;
				this.q[1 + qIndex] = 0 * this.q[0 + qIndex];
				this.q[2 + qIndex] = 0 * this.q[0 + qIndex];
				this.q[3 + qIndex] = 1 * this.q[0 + qIndex];
				qIndex += 4;
				xIndex += 2;
			}
		}
	},
	resetWave : function() {
		var xmid = .5 * (xmin + xmax);
		var ymid = .5 * (ymin + ymax);
		var dg = .2 * (xmax - xmin);
		var xIndex = 0;
		var qIndex = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.q[0 + qIndex] = .3 + Math.exp(-(Math.pow(this.x[0 + xIndex] - xmid, 2) + Math.pow(this.x[1 + xIndex] - ymid, 2))/(dg * dg));
				this.q[1 + qIndex] = 0 * this.q[0 + qIndex];
				this.q[2 + qIndex] = 0 * this.q[0 + qIndex];
				this.q[3 + qIndex] = 1 * this.q[0 + qIndex];
				xIndex += 2;
				qIndex += 4;
			}
		}
	},
	boundary : function() {
		this.boundaryMethod(this.nx, this.q);
	},
	step : function(dt) {
		
		//apply boundary conditions
		this.boundary();
	
		//solve
		this.advectMethod.advect.call(this, dt);
			
		//boundary again
		this.boundary();

		//compute pressure
		var qIndex = 0;
		var pIndex = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var rho = this.q[0 + qIndex];
				var u = this.q[1 + qIndex] / rho; 
				var v = this.q[2 + qIndex] / rho; 
				var energyTotal = this.q[3 + qIndex] / rho; 
				var energyKinematic = .5 * (u * u + v * v);
				var energyThermal = energyTotal - energyKinematic;
				this.pressure[pIndex] = (this.gamma - 1) * rho * energyThermal;
				++pIndex;
				qIndex += 4;
			}
		}
		
		//apply momentum diffusion = pressure
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var qIndex = 4 * (i + this.nx * j);
				for (var side = 0; side < 2; ++side) {
					var plusIndex = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					var minusIndex = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					this.q[1+side + qIndex] -= dt 
						* (this.pressure[plusIndex] - this.pressure[minusIndex]) 
						/ (this.x[side + 2 * plusIndex] - this.x[side + 2 * minusIndex]);
				}
			}
		}

		//apply work diffusion = momentum
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var qIndex = 4 * (i + this.nx * j);
				for (var side = 0; side < 2; ++side) {
					var qIndexR = 4 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]));
					var qIndexL = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
					//this is pulling the coordinate associated with the interface's direction
					//a more robust method would be to take both velocity components and dot them with the interface tangent
					var uR = this.q[1+side + qIndexR] / this.q[0 + qIndexR];
					var uL = this.q[1+side + qIndexL] / this.q[0 + qIndexL];
					var ipR = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					var ipL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					var ixR = side + 2 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]));
					var ixL = side + 2 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
					var dx = this.x[ixR] - this.x[ixL];
					this.q[3 + qIndex] -= dt * (this.pressure[ipR] * uR - this.pressure[ipL] * uL) / dx;
				}
			}
		}
	
		//last boundary update
		this.boundary();
	},
	update : function() {
		//get timestep
		var dt = this.advectMethod.initStep.call(this);

		//do the update
		this.step(dt);
	}
});


var Hydro = makeClass({
	init : function() {
		this.state = new HydroState();
	
		//geometry
		this.vertexPositions = new Float32Array(3*this.state.nx*this.state.nx);
		this.vertexStates = new Float32Array(3*this.state.nx*this.state.nx);
	
		//initialize geometry x and y coordinates once since they don't change
		var centerX = (xmax + xmin) / 2;
		var centerY = (ymax + ymin) / 2;
		var iv = 0;
		var ix = 0;
		var x = this.state.x;
		var nx = this.state.nx;
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				this.vertexPositions[iv] = x[ix] - centerX; ++iv; ++ix;
				this.vertexPositions[iv] = x[ix] - centerY; ++iv; ++ix;
				iv++;
			}
		}

	},
	update : function() {
		//todo adm or something
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();
	
		//update geometry z coordinate
		var q = this.state.q;
		var nx = this.state.nx;
		var iv = 2;
		var is = 0;
		var qIndex = 0;
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				this.vertexPositions[iv] = q[0 + qIndex] * 20.;
				iv += 3;
				this.vertexStates[is++] = q[0 + qIndex];
				this.vertexStates[is++] = Math.sqrt((q[1 + qIndex] * q[1 + qIndex] + q[2 + qIndex] * q[2 + qIndex]) / (q[0 + qIndex] * q[0 + qIndex]));
				this.vertexStates[is++] = q[3 + qIndex] / q[0 + qIndex];
				qIndex += 4;
			}
		}
	}
});

var hydro = new Hydro();

function update() {
	//iterate
	hydro.update();
	waveVtxBuf.updateData(hydro.vertexPositions);
	waveStateBuf.updateData(hydro.vertexStates);
	//draw
	GL.draw();
	requestAnimFrame(update);
}

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();
}

function buildSelect(id, key, map) {
	var select = $('#' + id);
	for (var k in map) {
		var option = $('<option>', {text : k});
		option.appendTo(select);
		if (hydro.state[key] == map[k]) {
			option.attr('selected', 'true');
		}
	}
	select.change(function() {
		hydro.state[key] = map[select.val()];
	});
}

$(document).ready(function(){
	panel = $('#panel');	

	$('#reset-sod').click(function(){ hydro.state.resetSod(); });
	$('#reset-wave').click(function(){ hydro.state.resetWave(); });

	buildSelect('boundary', 'boundaryMethod', boundaryMethods);
	buildSelect('flux-limiter', 'fluxMethod', fluxMethods);
	buildSelect('advect-method', 'advectMethod', advectMethods);

	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()
	
	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	gl.enable(gl.DEPTH_TEST);
	//GL.view.ortho = true;
	//GL.view.zNear = -1;
	//GL.view.zFar = 1;

	quat.multiply(GL.view.angle,
		[Math.sin(Math.PI/8), 0, 0, Math.cos(Math.PI/8)],
		[0,0,0,1]);

	var shader = new GL.ShaderProgram({
		vertexCodeID : 'water-vsh',
		vertexPrecision : 'best',
		fragmentCodeID : 'water-fsh',
		fragmentPrecision : 'best'
	});
	
	//make grid
	hydro.update();
	waveVtxBuf = new GL.ArrayBuffer({
		data : hydro.vertexPositions,
		usage : gl.DYNAMIC_DRAW
	});
	waveStateBuf = new GL.ArrayBuffer({
		data : hydro.vertexStates,
		usage : gl.DYNAMIC_DRAW
	});
	for (var j = 0; j < hydro.state.nx-1; ++j) {
		var indexes = [];
		for (var i = 0; i < hydro.state.nx; ++i) {
			indexes.push(i + j*hydro.state.nx);
			indexes.push(i + (j+1)*hydro.state.nx);
		}
		new GL.SceneObject({
			mode : gl.TRIANGLE_STRIP,
			indexes : new GL.ElementArrayBuffer({data:indexes}),
			attrs : {
				vertex : waveVtxBuf,
				state : waveStateBuf
			},
			shader : shader
		});
	}
	
	function refreshView() {
		GL.view.pos[0] = 0;
		GL.view.pos[1] = 0;
		GL.view.pos[2] = targetDistance;
		vec3.transformQuat(GL.view.pos, GL.view.pos, GL.view.angle);
	};
	var targetDistance = 2 * (xmin + xmax + ymin + ymax) / 4;
	var zoomFactor = .0003;	// upon mousewheel
	var dragging = false;
	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
			quat.multiply(GL.view.angle, GL.view.angle, tmpQ);
			quat.normalize(GL.view.angle, GL.view.angle);
			refreshView();
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-zoomFactor * zoomChange);
			targetDistance *= scale;
			refreshView();			
		}
	});
	refreshView();
	
	//start it off
	$(window).resize(onresize);
	onresize();
	update();
});
