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
function buildEigenstate(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, hTotal, gamma) {
	//calculate matrix & eigenvalues & vectors at interface from state at interface
	var speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * velocity * velocity));	
	//matrix, listed per column
	matrix[0][0] = 0;
	matrix[0][1] = (gamma - 3) / 2 * velocity * velocity;
	matrix[0][2] = velocity * ((gamma - 1) / 2 * velocity * velocity - hTotal);
	matrix[1][0] = 1;
	matrix[1][1] = (3 - gamma) * velocity;
	matrix[1][2] = hTotal - (gamma - 1) * velocity * velocity;
	matrix[2][0] = 0;
	matrix[2][1] = gamma - 1;
	matrix[2][2] = gamma * velocity;

	//eigenvalues: min, mid, max
	eigenvalues[0] = velocity - speedOfSound;
	eigenvalues[1] = velocity;
	eigenvalues[2] = velocity + speedOfSound;
	//min eigenvector
	eigenvectors[0][0] = 1;
	eigenvectors[0][1] = velocity - speedOfSound;
	eigenvectors[0][2] = hTotal - speedOfSound * velocity;
	//mid eigenvector
	eigenvectors[1][0] = 1;
	eigenvectors[1][1] = velocity;
	eigenvectors[1][2] = .5 * velocity * velocity;
	//max eigenvector
	eigenvectors[2][0] = 1;
	eigenvectors[2][1] = velocity + speedOfSound;
	eigenvectors[2][2] = hTotal + speedOfSound * velocity;
	//calculate eigenvector inverses ... 
	mat33invert(eigenvectorsInverse, eigenvectors);
}

var fluxMethods = {
	upwind : function(r) { return 0; },
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
			for (var iq = 0; iq < 4; ++iq) {
				//top
				q[iq + 4 * (i + nx * 0)] = q[iq + 4 * (i + nx * (nx-4))];
				q[iq + 4 * (i + nx * 1)] = q[iq + 4 * (i + nx * (nx-3))];
				q[iq + 4 * (i + nx * (nx-2))] = q[iq + 4 * (i + nx * 2)];
				q[iq + 4 * (i + nx * (nx-1))] = q[iq + 4 * (i + nx * 3)];
				//left	
				q[iq + 4 * (0 + nx * i)] = q[iq + 4 * (nx-4 + nx * i)];
				q[iq + 4 * (1 + nx * i)] = q[iq + 4 * (nx-3 + nx * i)];
				q[iq + 4 * (nx-2 + nx * i)] = q[iq + 4 * (2 + nx * i)];
				q[iq + 4 * (nx-1 + nx * i)] = q[iq + 4 * (3 + nx * i)];
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
			q[0 + 4 * (0 + nx * (i))] = q[0 + 4 * (3 + nx * (i))];
			q[0 + 4 * (1 + nx * (i))] = q[0 + 4 * (2 + nx * (i))];
			q[0 + 4 * (nx-2 + nx * (i))] = q[0 + 4 * (nx-3 + nx * (i))];
			q[0 + 4 * (nx-1 + nx * (i))] = q[0 + 4 * (nx-4 + nx * (i))];
			q[1 + 4 * (0 + nx * (i))] = -q[1 + 4 * (3 + nx * (i))];
			q[1 + 4 * (1 + nx * (i))] = -q[1 + 4 * (2 + nx * (i))];
			q[1 + 4 * (nx-2 + nx * (i))] = -q[1 + 4 * (nx-3 + nx * (i))];
			q[1 + 4 * (nx-1 + nx * (i))] = -q[1 + 4 * (nx-4 + nx * (i))];
			q[2 + 4 * (0 + nx * (i))] = -q[2 + 4 * (3 + nx * (i))];
			q[2 + 4 * (1 + nx * (i))] = -q[2 + 4 * (2 + nx * (i))];
			q[2 + 4 * (nx-2 + nx * (i))] = -q[2 + 4 * (nx-3 + nx * (i))];
			q[2 + 4 * (nx-1 + nx * (i))] = -q[2 + 4 * (nx-4 + nx * (i))];
			q[3 + 4 * (0 + nx * (i))] = q[3 + 4 * (3 + nx * (i))];
			q[3 + 4 * (1 + nx * (i))] = q[3 + 4 * (2 + nx * (i))];
			q[3 + 4 * (nx-2 + nx * (i))] = q[3 + 4 * (nx-3 + nx * (i))];
			q[3 + 4 * (nx-1 + nx * (i))] = q[3 + 4 * (nx-4 + nx * (i))];
		}
	},
	dirichlet : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var iq = 0; iq < 4; ++iq) {
				q[iq + 4 * (i + nx * (0))] = 0;
				q[iq + 4 * (1 + nx * (1))] = 0;
				q[iq + 4 * (i + nx * (nx-2))] = 0;
				q[iq + 4 * (i + nx * (nx-1))] = 0;
				q[iq + 4 * (0 + nx * (i))] = 0;
				q[iq + 4 * (1 + nx * (1))] = 0;
				q[iq + 4 * (nx-2 + nx * (i))] = 0;
				q[iq + 4 * (nx-1 + nx * (i))] = 0;
			}
		}
	},
	constant : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var iq = 0; iq < 4; ++iq) {
				q[iq + 4 * (i + nx * (0))] = q[iq + 4 * (i + nx * (1))] = q[iq + 4 * (i + nx * (2))];
				q[iq + 4 * (i + nx * (nx-1))] = q[iq + 4 * (i + nx * (nx-2))] = q[iq + 4 * (i + nx * (nx-3))];
				q[iq + 4 * (0 + nx * (i))] = q[iq + 4 * (1 + nx * (i))] = q[iq + 4 * (2 + nx * (i))];
				q[iq + 4 * (nx-1 + nx * (i))] = q[iq + 4 * (nx-2 + nx * (i))] = q[iq + 4 * (nx-3 + nx * (i))];
			}
		}
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {
			var mindum = undefined;
			for (var j = 0; j < this.nx; ++j) {
				for (var i = 0; i < this.nx; ++i) {
					var rho = this.q[0 + 4 * (i + this.nx * (j))];
					var u = this.q[1 + 4 * (i + this.nx * (j))] / rho;
					var v = this.q[2 + 4 * (i + this.nx * (j))] / rho; 
					var energyTotal = this.q[3 + 4 * (i + this.nx * (j))] / rho; 
					var energyKinematic = .5 * (u * u + v * v);
					var energyThermal = energyTotal - energyKinematic;
					var speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
					//just xi? yi too?
					var dum = (this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] - this.xi[0 + 2 * (i + (this.nx+1) * j)]) / (speedOfSound + Math.sqrt(u * u + v * v));
					if (mindum === undefined || dum < mindum) mindum = dum;
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
					//left boundary's du/dx (velocity gradient tangent to edge)
					this.ui[i][j][0] = .5 * (this.q[1 + 4 * (i + this.nx * (j))] / this.q[0 + 4 * (i + this.nx * (j))] + this.q[1 + 4 * (i-1 + this.nx * (j))] / this.q[0 + 4 * (i-1 + this.nx * (j))]);
					//top boundary's dv/dy
					this.ui[i][j][1] = .5 * (this.q[2 + 4 * (i + this.nx * (j))] / this.q[0 + 4 * (i + this.nx * (j))] + this.q[2 + 4 * (i + this.nx * (j-1))] / this.q[0 + 4 * (i + this.nx * (j-1))]);
				}
			}
			//boundary zero
			for (var j = 0; j < this.nghost; ++j) {
				for (var i = 0; i <= this.nx; ++i) {
					//left boundary, left and top sides, zero vector
					this.ui[j][i][0] = 0;
					this.ui[j][i][1] = 0;
					//right boundary, left and top sides, zero vector
					this.ui[this.nx-j][i][0] = 0;
					this.ui[this.nx-j][i][1] = 0;
					//top boundary, left and top sides, zero vector
					this.ui[i][j][0] = 0;
					this.ui[i][j][1] = 0;
					//bottom boundary, left and top sides, zero vector
					this.ui[i][this.nx-j][0] = 0;
					this.ui[i][this.nx-j][1] = 0;
				}
			}

			//compute flux and advect for each state vector
			for (var iq = 0; iq < 4; ++iq) {	//state
				//r_{i-1/2},{j-1/2} flux limiter
				for (var j = this.nghost; j < this.nx+this.nghost-3; ++j) {
					for (var i = this.nghost; i < this.nx+this.nghost-3; ++i) {
						
						//left boundary
						var dq = this.q[iq + 4 * (i + this.nx * (j))] - this.q[iq + 4 * (i-1 + this.nx * (j))];
						if (Math.abs(dq) > 0) {
							//left interface left velocity
							if (this.ui[i][j][0] >= 0) {
								this.r[i][j][0][iq] = (this.q[iq + 4 * (i-1 + this.nx * (j))] - this.q[iq + 4 * (i-2 + this.nx * (j))]) / dq;
							} else {
								this.r[i][j][0][iq] = (this.q[iq + 4 * (i+1 + this.nx * (j))] - this.q[iq + 4 * (i + this.nx * (j))]) / dq;
							}
						} else {
							this.r[i][j][0][iq] = 0;
						}
						
						//top boundary
						var dq = this.q[iq + 4 * (i + this.nx * (j))] - this.q[iq + 4 * (i + this.nx * (j-1))];
						if (Math.abs(dq) > 0) {
							//left interface left velocity
							if (this.ui[i][j][1] >= 0) {
								this.r[i][j][1][iq] = (this.q[iq + 4 * (i + this.nx * (j-1))] - this.q[iq + 4 * (i + this.nx * (j-2))]) / dq;
							} else {
								this.r[i][j][1][iq] = (this.q[iq + 4 * (i + this.nx * (j+1))] - this.q[iq + 4 * (i + this.nx * (j))]) / dq;
							}
						} else {
							this.r[i][j][1][iq] = 0;
						}
					}
				}
			
				//now for ghost boundaries
				for (var j = 0; j < this.nghost; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							this.r[i][j][side][iq] = 0;
							this.r[i][this.nx-j][side][iq] = 0;
							this.r[j][i][side][iq] = 0;
							this.r[this.nx-j][i][side][iq] = 0;
						}
					}
				}

				//construct flux:
				for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						
						//left
						//flux limiter
						var phi = this.fluxMethod(this.r[i][j][0][iq]);
						if (this.ui[i][j][0] >= 0) {
							this.flux[i][j][0][iq] = this.ui[i][j][0] * this.q[iq + 4 * (i-1 + this.nx * (j))];
						} else {
							this.flux[i][j][0][iq] = this.ui[i][j][0] * this.q[iq + 4 * (i + this.nx * (j))];
						}
						var delta = phi * (this.q[iq + 4 * (i + this.nx * (j))] - this.q[iq + 4 * (i-1 + this.nx * (j))]);
						var dx = this.x[0 + 2 * (i + this.nx * j)] - this.x[0 + 2 * (i-1 + this.nx * j)];
						this.flux[i][j][0][iq] += delta * .5 * Math.abs(this.ui[i][j][0]) * (1 - Math.abs(this.ui[i][j][0] * dt / dx));
					
						//top
						//flux limiter
						var phi = this.fluxMethod(this.r[i][j][1][iq]);
						if (this.ui[i][j][1] >= 0) {
							this.flux[i][j][1][iq] = this.ui[i][j][1] * this.q[iq + 4 * (i + this.nx * (j-1))];
						} else {
							this.flux[i][j][1][iq] = this.ui[i][j][1] * this.q[iq + 4 * (i + this.nx * (j))];
						}
						var delta = phi * (this.q[iq + 4 * (i + this.nx * (j))] - this.q[iq + 4 * (i + this.nx * (j-1))]);
						var dx = this.x[1 + 2 * (i + this.nx * j)] - this.x[1 + 2 * (i + this.nx * (j-1))];
						this.flux[i][j][1][iq] += delta * .5 * Math.abs(this.ui[i][j][1]) * (1 - Math.abs(this.ui[i][j][1] * dt / dx));
					}
				}
				//now for ghost boundaries
				for (var j = 0; j < this.nghost-1; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							this.flux[i][j][side][iq] = 0;
							this.flux[i][this.nx-j][side][iq] = 0;
							this.flux[j][i][side][iq] = 0;
							this.flux[this.nx-j][i][side][iq] = 0;
						}
					}
				}

				//update cells
				for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
						this.q[iq + 4 * (i + this.nx * (j))] -= dt * (this.flux[i+1][j][0][iq] - this.flux[i][j][0][iq]) / (this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] - this.xi[0 + 2 * (i + (this.nx+1) * j)]);
						this.q[iq + 4 * (i + this.nx * (j))] -= dt * (this.flux[i][j+1][1][iq] - this.flux[i][j][1][iq]) / (this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] - this.xi[1 + 2 * (i + (this.nx+1) * j)]);
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

		//q_i,j,state: state vector, stored as q[state + 4 * (j + this.nx * (i))]
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
		
		//p_i: pressure
		this.pressure = [];//new Float32Array(this.nx);
		for (var i = 0; i < this.nx; ++i) {
			this.pressure[i] = [];
			for (var j = 0; j < this.nx; ++j) {
				this.pressure[i][j] = 0;
			}
		}

		//used for Burgers
		//r_{i-1/2},{j-1/2},side,state	
		this.r = [];
		for (var i = 0; i < this.nx+1; ++i) {
			this.r[i] = [];
			for (var j = 0; j < this.nx+1; ++j) {
				this.r[i][j] = [];
				for (var side = 0; side < 2; ++side) {
					this.r[i][j][side] = [0,0,0,0];
				}
			}
		}
		
		//f_{i-1/2},{j-1/2}: cell flux
		this.flux = [];
		for (var i = 0; i < this.nx+1; ++i) {
			this.flux[i] = [];
			for (var j = 0; j < this.nx+1; ++j) {
				this.flux[i][j] = [
					[0,0,0,0],
					[0,0,0,0]
				];
			}
		}
		
		//only used with Burger's eqn advection code
		//u_{i-1/2}: interface velocity
		this.ui = [];//new Float32Array(this.nx+1);
		for (var i = 0; i < this.nx+1; ++i) {
			this.ui[i] = [];
			for (var j = 0; j < this.nx+1; ++j) {
				this.ui[i][j] = [];
				for (var side = 0; side < 2; ++side) {
					this.ui[i][j][side] = 0;
				}
			}
		}

		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethod = boundaryMethods.mirror;
		this.fluxMethod = fluxMethods.upwind;
		this.advectMethod = advectMethods.Burgers;
	},
	resetSod : function() {
		var ix = 0;
		var iq = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.q[0 + iq] = (this.x[0 + ix] < 30 && this.x[1 + ix] < 30) ? 1 : .1;
				this.q[1 + iq] = 0 * this.q[0 + iq];
				this.q[2 + iq] = 0 * this.q[0 + iq];
				this.q[3 + iq] = 1 * this.q[0 + iq];
				iq += 4;
				ix += 2;
			}
		}
	},
	resetWave : function() {
		var xmid = .5 * (xmin + xmax);
		var ymid = .5 * (ymin + ymax);
		var dg = .2 * (xmax - xmin);
		var ix = 0;
		var iq = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				this.q[0 + iq] = .3 + Math.exp(-(Math.pow(this.x[0 + ix] - xmid, 2) + Math.pow(this.x[1 + ix] - ymid, 2))/(dg * dg));
				this.q[1 + iq] = 0 * this.q[0 + iq];
				this.q[2 + iq] = 0 * this.q[0 + iq];
				this.q[3 + iq] = 1 * this.q[0 + iq];
				ix += 2;
				iq += 4;
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
		var iq = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var rho = this.q[0 + iq];
				var u = this.q[1 + iq] / rho; 
				var v = this.q[2 + iq] / rho; 
				var energyTotal = this.q[3 + iq] / rho; 
				var energyKinematic = .5 * (u * u + v * v);
				var energyThermal = energyTotal - energyKinematic;
				this.pressure[i][j] = (this.gamma - 1) * rho * energyThermal;
				iq += 4;
			}
		}
		
		//apply momentum diffusion = pressure
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				this.q[1 + 4 * (i + this.nx * (j))] -= dt * (this.pressure[i+1][j] - this.pressure[i-1][j]) / (this.x[0 + 2 * (i+1 + this.nx * j)] - this.x[0 + 2 * (i-1 + this.nx * j)]);
				this.q[2 + 4 * (i + this.nx * (j))] -= dt * (this.pressure[i][j+1] - this.pressure[i][j-1]) / (this.x[1 + 2 * (i + this.nx * (j+1))] - this.x[1 + 2 * (i + this.nx * (j-1))]);
			}
		}

		//apply work diffusion = momentum
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				//perhaps ui isn't the best use for velocit here?
				//perhaps we could derive something from the state variables?
				var u_inext = this.q[1 + 4 * (i+1 + this.nx * (j))] / this.q[0 + 4 * (i+1 + this.nx * (j))];
				var u_iprev = this.q[1 + 4 * (i-1 + this.nx * (j))] / this.q[0 + 4 * (i-1 + this.nx * (j))];
				var v_inext = this.q[2 + 4 * (i + this.nx * (j+1))] / this.q[0 + 4 * (i + this.nx * (j+1))];
				var v_iprev = this.q[2 + 4 * (i + this.nx * (j-1))] / this.q[0 + 4 * (i + this.nx * (j-1))];
				this.q[3 + 4 * (i + this.nx * (j))] -= dt * (this.pressure[i+1][j] * u_inext - this.pressure[i-1][j] * u_iprev) / (this.x[0 + 2 * (i+1 + this.nx * j)] - this.x[0 + 2 * (i-1 + this.nx * j)]);
				this.q[3 + 4 * (i + this.nx * (j))] -= dt * (this.pressure[i][j+1] * v_inext - this.pressure[i][j-1] * v_iprev) / (this.x[1 + 2 * (i + this.nx * (j+1))] - this.x[1 + 2 * (i + this.nx * (j-1))]);
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
		var iq = 0;
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				this.vertexPositions[iv] = q[0 + iq] * 20.;
				iv += 3;
				this.vertexStates[is++] = q[0 + iq];
				this.vertexStates[is++] = Math.sqrt((q[1 + iq] * q[1 + iq] + q[2 + iq] * q[2 + iq]) / (q[0 + iq] * q[0 + iq]));
				this.vertexStates[is++] = q[3 + iq] / q[0 + iq];
				iq += 4;
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
	onresize();
	update();
});
