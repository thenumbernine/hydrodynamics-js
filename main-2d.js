/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured

sources:
http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/
http://www.cfdbooks.com/cfdcodes.html
"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro
http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf also does not
*/

var panel;
var canvas;
var waveVtxBuf, waveStateBuf;
var xmin = -.5;
var xmax = .5; 
var ymin = -.5;
var ymax = .5;
var useNoise = true;
var mouse;

//interface directions
var dirs = [[1,0], [0,1]];

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

	//I'm going with http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf

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
	var identCheck = [];
	var identBad = false;
	for (var i = 0; i < 4; ++i) {
		for (var j = 0; j < 4; ++j) {
			var s = 0;
			var d = 0;
			for (var k = 0; k < 4; ++k) {
				/** /
				s += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)] * eigenvalues[k + 4 * offset];
				identCheck += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)];
				/**/
				s += eigenvectors[i + 4 * (k + 4 * offset)] * eigenvalues[k + 4 * offset] * eigenvectorsInverse[k + 4 * (j + 4 * offset)];
				d += eigenvectors[i + 4 * (k + 4 * offset)] * eigenvectorsInverse[k + 4 * (j + 4 * offset)];
				/**/
			}
			matrix[i + 4 * (j + 4 * offset)] = s;
			identCheck[i + 4 * j] = d;
			var epsilon = 1e-5;
			if (Math.abs(d - (i == j ? 1 : 0)) > epsilon) identBad = true;
		}
	}
	if (identBad) {
		console.log('bad eigen basis', identCheck);
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
	none : function() {},
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
					var energyKinetic = .5 * (u * u + v * v);
					var energyThermal = energyTotal - energyKinetic;
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
					var uiIndex = 2 * (i + (this.nx+1) * j);
					for (var side = 0; side < 2; ++side) {
						var qIndexR = 4 * (i + this.nx * j);
						var qIndexL = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
						this.ui[side + uiIndex] = .5 * (
							this.q[1+side + qIndexR] / this.q[0 + qIndexR] + 
							this.q[1+side + qIndexL] / this.q[0 + qIndexL]);
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


				var dxi = [];
				//construct flux:
				for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						var dx = this.x[0 + 2 * (i + this.nx * j)] 
							- this.x[0 + 2 * (i-dirs[0][0] + this.nx * (j-dirs[0][1]))];
						var dy = this.x[1 + 2 * (i + this.nx * j)] 
							- this.x[1 + 2 * (i-dirs[1][0] + this.nx * (j-dirs[1][1]))];
						dxi[0] = dx;
						dxi[1] = dy;
						var volume = dx * dy;
						//flux calculation 
						for (var side = 0; side < 2; ++side) {
							
							var rIndex = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							//apply limiter
							var phi = this.fluxMethod(this.r[rIndex]);
							var uiIndex = side + 2 * (i + (this.nx+1) * j);
							var fluxIndex = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							var qIndexL = state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
							var qIndexR = state + 4 * (i + this.nx * j);
							if (this.ui[uiIndex] >= 0) {
								this.flux[fluxIndex] = this.ui[uiIndex] * this.q[qIndexL];
							} else {
								this.flux[fluxIndex] = this.ui[uiIndex] * this.q[qIndexR];
							}
							var delta = phi * (this.q[qIndexR] - this.q[qIndexL]);
							this.flux[fluxIndex] += delta * .5
								* Math.abs(this.ui[uiIndex])
								* (1 - Math.abs(this.ui[uiIndex] * dt / dxi[side]));//* volume / (dxi[side] * dxi[side])));
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
				var dxi = [];
				for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
						
						var xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
						var xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
						var dx = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						var xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
						var xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
						var dy = this.xi[xiIndexR] - this.xi[xiIndexL];
					
						var volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						
						for (var side = 0; side < 2; ++side) {
							var ifluxR = state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
							var ifluxL = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							var df = this.flux[ifluxR] - this.flux[ifluxL];
							this.q[state + 4 * (i + this.nx * j)] -= dt * df / dxi[side];//* volume / (dxi[side] * dxi[side]);
						}
					}
				}
			}
		}
	},
	Riemann : {
		initStep : function() {
			var mindum = undefined;
			for (var j = 1; j < this.nx; ++j) {
				for (var i = 1; i < this.nx; ++i) {
					for (var side = 0; side < 2; ++side) {
						var qIndexL = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
						var densityL = this.q[0 + qIndexL];
						var velocityXL = this.q[1 + qIndexL] / densityL;
						var velocityYL = this.q[2 + qIndexL] / densityL;
						var energyTotalL = this.q[3 + qIndexL] / densityL;
						var energyKineticL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
						var energyThermalL = energyTotalL - energyKineticL;
						var pressureL = (this.gamma - 1) * densityL * energyThermalL;
						var speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
						var hTotalL = energyTotalL + pressureL / densityL;
						var roeWeightL = Math.sqrt(densityL);
						
						var qIndexR = 4 * (i + this.nx * j);
						var densityR = this.q[0 + qIndexR];
						var velocityXR = this.q[1 + qIndexR] / densityR;
						var velocityYR = this.q[2 + qIndexR] / densityR;
						var energyTotalR = this.q[3 + qIndexR] / densityR;
						var energyKineticR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
						var energyThermalR = energyTotalR - energyKineticR;
						var pressureR = (this.gamma - 1) * densityR * energyThermalR;
						var speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
						var hTotalR = energyTotalR + pressureR / densityR;
						var roeWeightR = Math.sqrt(densityR);

						var denom = roeWeightL + roeWeightR;
						var velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
						var velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
						var hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
						
						buildEigenstate(
							 //index into interface element.  
							 //from there you'll have to scale by cell size.  
							 //Thus manually recreating the automatic storage of C structures. 
							 //JavaScript, why can't you be more like LuaJIT? 
							 side + 2 * (i + (this.nx+1) * j),	
							 this.interfaceMatrix,	//dim^2 = 16
							 this.interfaceEigenvalues,	//dim = 4
							 this.interfaceEigenvectors,	//dim^2 = 16
							 this.interfaceEigenvectorsInverse,	//dim^2 = 16
							 velocityX, velocityY, hTotal, this.gamma,
							 dirs[side][0], dirs[side][1]);

						var maxLambda = Math.max(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[1+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[2+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[3+4*(side+2*(i+(this.nx+1)*j))]);
						var minLambda = Math.min(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[1+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[2+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[3+4*(side+2*(i+(this.nx+1)*j))]);
						var dx = this.xi[side + 2 * (i+dirs[side][0] + (this.nx+1) * (j+dirs[side][1]))] 
							- this.xi[side + 2 * (i + (this.nx+1) * j)];
						var dum = dx / (maxLambda - minLambda);
						if (mindum === undefined || dum < mindum) mindum = dum;
					}
				}
			}
			return this.cfl * mindum;
	
		},
		advect : function(dt) {
			for (var j = 1; j < this.nx; ++j) {
				for (var i = 1; i < this.nx; ++i) {
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							//find state change across interface in the basis of the eigenspace at the interface
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
									//reproject into interface eigenspace
								sum += this.interfaceEigenvectorsInverse[state + 4 * (k + 4 * (side + 2 * (i + (this.nx+1) * j)))]
									//flux difference
									* (this.q[k + 4 * (i + this.nx * j)] 
										- this.q[k + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))])
							}
							this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = sum;
						}
					}
				}
			}
			
			//boundary zero
			for (var j = 0; j < this.nghost-1; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					for (var state = 0; state < 4; ++state) {
						//left boundary, left and top sides, zero vector
						this.interfaceDeltaQTilde[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
						//right boundary, left and top sides, zero vector
						this.interfaceDeltaQTilde[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						//top boundary, left and top sides, zero vector
						this.interfaceDeltaQTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
						//bottom boundary, left and top sides, zero vector
						this.interfaceDeltaQTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
					}
				}
			}

			for (var j = this.nghost; j < this.nx + this.nghost - 3; ++j) {
				for (var i = this.nghost; i < this.nx + this.nghost - 3; ++i) {
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							var interfaceDeltaQTilde = this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))];
							if (Math.abs(interfaceDeltaQTilde) > 0) {
								if (this.interfaceEigenvalues[state + 4 * (side + 2 * (i + (this.nx+1) * j))] > 0) {
									this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i - dirs[side][0] + (this.nx+1) * (j - dirs[side][1])))]
										/ interfaceDeltaQTilde;
								} else {
									this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))]
										/ interfaceDeltaQTilde;
								}
							} else {
								this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							}
						}
					}
				}
			}
	
			//..and keep the boundary rTilde's zero	
			for (var j = 0; j < this.nghost; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					for (var state = 0; state < 4; ++state) {
						//left
						this.rTilde[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
						//right
						this.rTilde[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						//bottom
						this.rTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
						//top
						this.rTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
					}
				}
			}
		
			var fluxAvg = [];	//4
			var fluxTilde = [];	//4
			var dxi = [];
			//transform cell q's into cell qTilde's (eigenspace)
			for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
				for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
					var dx = this.xi[0 + 2 * (i + (this.nx+1) * j)] 
						- this.xi[0 + 2 * (i-dirs[0][0] + (this.nx+1) * (j-dirs[0][1]))];
					var dy = this.xi[1 + 2 * (i + (this.nx+1) * j)] 
						- this.xi[1 + 2 * (i-dirs[1][0] + (this.nx+1) * (j-dirs[1][1]))];
					var volume = dx * dy;
					dxi[0] = dx;
					dxi[1] = dy;
					for (var side = 0; side < 2; ++side) {
						
						//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
						//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
						for (var state = 0; state < 4; ++state) {
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
								sum += this.interfaceMatrix[state + 4 * (k + 4 * (side + 2 * (i + (this.nx+1) * j)))]
									* (this.q[k + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))]
										+ this.q[k + 4 * (i + this.nx * j)]);
							}
							fluxAvg[state] = .5 * sum;
						}

						//calculate flux
						for (var state = 0; state < 4; ++state) {
							var theta = 0;
							if (this.interfaceEigenvalues[state + 4 * (side + 2 * (i + (this.nx+1) * j))] >= 0) {
								theta = 1;
							} else {
								theta = -1;
							}
						
							var phi = this.fluxMethod(this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))]);

							var epsilon = this.interfaceEigenvalues[state + 4 * (side + 2 * (i + (this.nx+1) * j))] * dt / dxi[side];//* volume / (dxi[side] * dxi[side]); 

							var deltaFluxTilde = this.interfaceEigenvalues[state + 4 * (side + 2 * (i + (this.nx+1) * j))]
								* this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))];

							fluxTilde[state] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
						}
					
						//reproject fluxTilde back into q
						for (var state = 0; state < 4; ++state) {
							var sum = 0;
							for (var k = 0; k < 4; ++k) {
								sum += fluxTilde[k]
									* this.interfaceEigenvectors[state + 4 * (k + 4 * (side + 2 * (i + (this.nx+1) * j)))];
							}
							this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = fluxAvg[state] + sum;
						}
					}
				}
			}
		
			//zero boundary flux
			//..and keep the boundary r's zero	
			for (var j = 0; j < this.nghost-1; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					for (var state = 0; state < 4; ++state) {
						//left
						this.flux[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
						//right
						this.flux[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						//bottom
						this.flux[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
						//top
						this.flux[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
					}
				}
			}

			//update cells
			for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
					var xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
					var xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
					var dx = this.xi[xiIndexR] - this.xi[xiIndexL];
					
					var xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
					var xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
					var dy = this.xi[xiIndexR] - this.xi[xiIndexL];
					
					var volume = dx * dy;
					dxi[0] = dx;
					dxi[1] = dy;
					
					for (var side = 0; side < 2; ++side) {
						for (var state = 0; state < 4; ++state) {
							
							var ifluxR = state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
							var ifluxL = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							var df = this.flux[ifluxR] - this.flux[ifluxL];
							this.q[state + 4 * (i + this.nx * j)] -= dt * df / dxi[side];//* volume / (dxi[side] * dxi[side]);
						}
					}
				}
			}
		}
	}
};

var HydroState = makeClass({ 
	init : function(args) {
		this.nx = args.size;
		this.cfl =.5;
		this.gamma = args.gamma;
		
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
		this.advectMethod = advectMethods.Riemann;
	},
	resetSod : function() {
		var xIndex = 0;
		var qIndex = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var x = this.x[0 + xIndex];
				var y = this.x[1 + xIndex];
				var rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
				var u = 0;
				var v = 0;
				if (useNoise) {
					u += (Math.random() - .5) * 2 * .01;
					v += (Math.random() - .5) * 2 * .01;
				}
				var energyKinetic = .5 * (u * u + v * v);
				var energyThermal = 1;
				var energyTotal = energyKinetic + energyThermal;
				this.q[0 + qIndex] = rho;
				this.q[1 + qIndex] = rho * u; 
				this.q[2 + qIndex] = rho * v; 
				this.q[3 + qIndex] = rho * energyTotal; 
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
				var x = this.x[0 + xIndex];
				var y = this.x[1 + xIndex];
				var dx = x - xmid;
				var dy = y - ymid;
				var rho = 3 * Math.exp(-(dx * dx + dy * dy) / (dg * dg)) + .1;
				var u = 0; 
				var v = 0;
				if (useNoise) {
					u += (Math.random() - .5) * 2 * .01;
					v += (Math.random() - .5) * 2 * .01;
				}
				var energyKinetic = .5 * (u * u + v * v);
				var energyThermal = 1;
				var energyTotal = energyKinetic + energyThermal;
				this.q[0 + qIndex] = rho;
				this.q[1 + qIndex] = rho * u;
				this.q[2 + qIndex] = rho * v;
				this.q[3 + qIndex] = rho * energyTotal; 
				xIndex += 2;
				qIndex += 4;
			}
		}
	},
	//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
	resetKelvinHemholtz : function() {
		var xmid = .5 * (xmin + xmax);
		var xIndex = 0;
		var qIndex = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var x = this.x[0 + xIndex];
				var y = this.x[1 + xIndex];
				var yInTheMiddle = y > (.75 * ymin + .25 * ymax) && y < (.25 * ymin + .75 * ymax);
				var rho = yInTheMiddle ? 2 : 1;
				var u = yInTheMiddle ? .5 : -.5;
				var v = 0;
				if (useNoise) {
					u += (Math.random() - .5) * 2 * .01;
					v += (Math.random() - .5) * 2 * .01;
				}
				//P = (gamma - 1) rho (eTotal - eKinetic)
				//eTotal = P / ((gamma - 1) rho) + eKinetic
				var pressure = 2.5;
				var energyKinetic = .5 * (u * u + v * v);
				var energyTotal = pressure / ((this.gamma - 1) * rho) + energyKinetic; 
				this.q[0 + qIndex] = rho;
				this.q[1 + qIndex] = rho * u; 
				this.q[2 + qIndex] = rho * v;
				this.q[3 + qIndex] = rho * energyTotal; 
				xIndex += 2;
				qIndex += 4;
			}
		}
		//TODO make it periodic on the left/right borders and reflecting on the top/bottom borders
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
				var energyKinetic = .5 * (u * u + v * v);
				var energyThermal = energyTotal - energyKinetic;
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
		var dxi = [];
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var qIndex = 4 * (i + this.nx * j);
				for (var side = 0; side < 2; ++side) {
					var plusIndex = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					var minusIndex = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					//this is pulling the coordinate associated with the interface's direction
					//a more robust method would be to take both velocity components and dot them with the interface tangent
					var uR = this.q[1+side + 4 * plusIndex] / this.q[0 + 4 * plusIndex];
					var uL = this.q[1+side + 4 * minusIndex] / this.q[0 + 4 * minusIndex];
					this.q[3 + qIndex] -= dt * (this.pressure[plusIndex] * uR - this.pressure[minusIndex] * uL) 
								/ (this.x[side + 2 * plusIndex] - this.x[side + 2 * minusIndex]);
				}
			}
		}
	
		//last boundary update
		this.boundary();
	},
	update : function() {
		//get timestep
		var dt = this.advectMethod.initStep.call(this);
window.lastDT = dt;
		//do the update
		this.step(dt);
	}
});


var Hydro = makeClass({
	init : function() {
		this.state = new HydroState({
			size : 200, 
			gamma : 7/5
		});
	
		//geometry
		this.vertexPositions = new Float32Array(2*this.state.nx*this.state.nx);
		this.vertexStates = new Float32Array(this.state.nx*this.state.nx);
	
		//initialize geometry x and y coordinates once since they don't change
		var centerX = (xmax + xmin) / 2;
		var centerY = (ymax + ymin) / 2;
		var vIndex = 0;
		var xIndex = 0;
		var x = this.state.x;
		var nx = this.state.nx;
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				this.vertexPositions[vIndex] = x[xIndex] - centerX; ++vIndex; ++xIndex;
				this.vertexPositions[vIndex] = x[xIndex] - centerY; ++vIndex; ++xIndex;
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
		var pressure = this.state.pressure;
		var nx = this.state.nx;
		var e = 0;
	
		//density
		var dataMin = q[0];
		var dataMax = q[0];
		//pressure
		//var dataMin = pressure[0];
		//var dataMax = pressure[0] + 1e-9;
		
		var lastDataRange = this.lastDataMax - this.lastDataMin;	
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				//density
				var s = q[0 + 4 * e];
				//pressure
				//var s = pressure[e];
				
				if (s < dataMin) dataMin = s;
				if (s > dataMax) dataMax = s;
				
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		if (this.updateLastDataRange) {
			this.lastDataMin = dataMin;
			this.lastDataMax = dataMax;
			var thisTick = Math.floor(Date.now() / 1000);
			if (this.lastDataRangeUpdateTime != thisTick) {
				this.lastDataRangeUpdateTime = thisTick;
				$('#dataRangeFixedMin').val(this.lastDataMin);
				$('#dataRangeFixedMax').val(this.lastDataMax);
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

var sceneObjects = [];
var colorSchemes = {};

$(document).ready(function(){
	panel = $('#panel');	

	$('#reset-sod').click(function(){ hydro.state.resetSod(); });
	$('#reset-wave').click(function(){ hydro.state.resetWave(); });
	$('#reset-kelvin-hemholtz').click(function(){ hydro.state.resetKelvinHemholtz(); });

	$('#use-noise').change(function() {
		useNoise = $(this).is(':checked');
	});

	buildSelect('boundary', 'boundaryMethod', boundaryMethods);
	buildSelect('flux-limiter', 'fluxMethod', fluxMethods);
	buildSelect('advect-method', 'advectMethod', advectMethods);

	hydro.lastDataMin = Number($('#dataRangeFixedMin').val());
	hydro.lastDataMax = Number($('#dataRangeFixedMax').val());
	hydro.updateLastDataRange = true;
	$('#dataRangeScaleNormalized').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = true;
	});
	$('#dataRangeScaleFixed').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = false;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});
	$('#dataRangeFixedMin').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
	});
	$('#dataRangeFixedMax').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});

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
	GL.view.ortho = true;
	GL.view.zNear = -1;
	GL.view.zFar = 1;
	GL.view.fovY = 125 / 200 * (xmax - xmin);

	colorSchemes.Heat = new GL.GradientTexture({
		width:256, 
		colors:[
			[0,0,.5],
			[0,0,1],
			[1,1,0],
			[1,0,0],
		],
		dontRepeat : true
	});

	var isobarSize = 16;
	var isobarData = new Uint8Array(isobarSize);
	for (var i = 1; i < isobarSize; i += 2) {
		isobarData[i] = 255;
	}
	colorSchemes['B&W'] = new GL.Texture2D({
		width : isobarSize,
		height : 1,
		format : gl.LUMINANCE,
		internalFormat : gl.LUMINANCE,
		data : isobarData,
		minFilter : gl.LINEAR,
		magFilter : gl.NEAREST,
		generateMipmap : true,
		wrap : {s : gl.REPEAT, t : gl.CLAMP_TO_EDGE }
	});

	colorSchemes.Heat.bind();

	for (var k in colorSchemes) {
		(function(){
			var v = colorSchemes[k];
			$('#color-scheme').append($('<option>', {
				value : k,
				text : k
			}));
		})();
	}
	$('#color-scheme').change(function() {
		var k = $(this).val();
		var v = colorSchemes[k];
		$.each(sceneObjects, function(k, sceneObject) {
			gl.bindTexture(gl.TEXTURE_2D, v.obj);
		});
	});

	var shader = new GL.ShaderProgram({
		vertexCodeID : 'water-vsh',
		vertexPrecision : 'best',
		fragmentCodeID : 'water-fsh',
		fragmentPrecision : 'best'
	});
	
	//make grid
	hydro.update();
	waveVtxBuf = new GL.ArrayBuffer({
		dim : 2,
		data : hydro.vertexPositions,
		usage : gl.DYNAMIC_DRAW
	});
	waveStateBuf = new GL.ArrayBuffer({
		dim : 1,
		data : hydro.vertexStates,
		usage : gl.DYNAMIC_DRAW
	});
	for (var j = 0; j < hydro.state.nx-1; ++j) {
		var indexes = [];
		for (var i = 0; i < hydro.state.nx; ++i) {
			indexes.push(i + j*hydro.state.nx);
			indexes.push(i + (j+1)*hydro.state.nx);
		}
		sceneObjects.push(new GL.SceneObject({
			mode : gl.TRIANGLE_STRIP,
			indexes : new GL.ElementArrayBuffer({data:indexes}),
			attrs : {
				vertex : waveVtxBuf,
				state : waveStateBuf
			},
			uniforms : {
				tex : 0
			},
			shader : shader
		}));
	}
	
	var zoomFactor = .0003;	// upon mousewheel
	var dragging = false;
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			var aspectRatio = canvas.width / canvas.height;
			GL.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * GL.view.fovY);
			GL.view.pos[1] += dy / canvas.height * 2 * GL.view.fovY;
			GL.updateProjection();
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-zoomFactor * zoomChange);
			GL.view.fovY *= scale 
			GL.updateProjection();
		}
	});
	
	//start it off
	$(window).resize(onresize);
	onresize();
	update();
});
