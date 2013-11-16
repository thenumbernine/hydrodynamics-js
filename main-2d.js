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
	//calculate matrix & eigenvalues & vectors at interface from state at interface
	var speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)));
	var velocityN = velocityX * normalX + velocityY * normalY;
	var tangentX = -normalY;
	var tangentY = normalX;
	var velocityT = velcoityX * tangentX + velocityY * tangentY;
	
	//eigenvalues: min, mid, max
	eigenvalues[0] = velocityN - speedOfSound;
	eigenvalues[1] = velocityN;
	eigenvalues[2] = velocityN;
	eigenvalues[3] = velocityN + speedOfSound;
	
	//min eigenvector
	eigenvectors[0][0] = 1;
	eigenvectors[0][1] = velocityX - speedOfSound * normalX;
	eigenvectors[0][2] = velocityY - speedOfSound * normalY;
	eigenvectors[0][3] = hTotal - speedOfSound * velocityN;
	//mid eigenvector
	eigenvectors[1][0] = 0;
	eigenvectors[1][1] = speedOfSound * tangentX;
	eigenvectors[1][2] = speedOfSound * tangentY;
	eigenvectors[1][3] = speedOfSound * velocityT;
	//mid eigenvector
	eigenvectors[2][0] = 1;
	eigenvectors[2][1] = velocityX;
	eigenvectors[2][2] = velocityY;
	eigenvectors[2][3] = .5 * (velocityX * velocityX + velocityY * velocityY);
	//max eigenvector
	eigenvectors[3][0] = 1;
	eigenvectors[3][1] = velocityX + speedOfSound * normalX;
	eigenvectors[3][2] = velocityY + speedOfSound * normalY;
	eigenvectors[3][3] = hTotal + speedOfSound * velocityN;
	
	//calculate eigenvector inverses ... 
	mat44invert(eigenvectorsInverse, eigenvectors);
	
	//calculate matrix
	matrix[0][0] = 0;
	matrix[0][1] = (gamma - 3) / 2 * velocityX * velocityX + (gamma - 1) / 2 * velocityY * velocityY;
	matrix[0][2] = -velocityX * velocityY;
	matrix[0][3] = (gamma - 1) * q * q * u - gamma * E * u / rho;//velocity * ((gamma - 1) / 2 * velocity * velocity - hTotal);
	matrix[1][0] = 1;
	matrix[1][1] = (3 - gamma) * velocityX;
	matrix[1][2] = velocityY;
	matrix[1][3] = (gamma - 1) / 2 * (3 * velocityX * velocityX + velocityY * velocityY);//hTotal - (gamma - 1) * velocity * velocity;
	matrix[2][0] = 0;
	matrix[2][1] = (1 - gamma) * velocityY;
	matrix[2][2] = velocityX;
	matrix[2][3] = (1 - gamma) * velocityX * velocityY;
	matrix[2][0] = 0;
	matrix[2][1] = gamma - 1;
	matrix[2][2] = 0;
	matrix[2][3] = gamma * velocityX;


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
			for (var iq = 0; iq < 4; ++iq) {
				q[iq + 4 * (i + nx * (0))] = 0;
				q[iq + 4 * (1 + nx * (1))] = 0;
				q[iq + 4 * (i + nx * (nx-2))] = 0;
				q[iq + 4 * (i + nx * (nx-1))] = 0;
				q[iq + 4 * (0 + nx * i)] = 0;
				q[iq + 4 * (1 + nx * (1))] = 0;
				q[iq + 4 * (nx-2 + nx * i)] = 0;
				q[iq + 4 * (nx-1 + nx * i)] = 0;
			}
		}
	},
	constant : function(nx,q) {
		for (var i = 0; i < nx; ++i) {
			for (var iq = 0; iq < 4; ++iq) {
				q[iq + 4 * (i + nx * (0))] = q[iq + 4 * (i + nx * (1))] = q[iq + 4 * (i + nx * (2))];
				q[iq + 4 * (i + nx * (nx-1))] = q[iq + 4 * (i + nx * (nx-2))] = q[iq + 4 * (i + nx * (nx-3))];
				q[iq + 4 * (0 + nx * i)] = q[iq + 4 * (1 + nx * i)] = q[iq + 4 * (2 + nx * i)];
				q[iq + 4 * (nx-1 + nx * i)] = q[iq + 4 * (nx-2 + nx * i)] = q[iq + 4 * (nx-3 + nx * i)];
			}
		}
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {
			var mindum = undefined;
			var iq = 0;
			for (var j = 0; j < this.nx; ++j) {
				for (var i = 0; i < this.nx; ++i) {
					var rho = this.q[0 + iq];
					var u = this.q[1 + iq] / rho;
					var v = this.q[2 + iq] / rho; 
					var energyTotal = this.q[3 + iq] / rho; 
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
					iq += 4;
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
			for (var iq = 0; iq < 4; ++iq) {	//state
				//r_{i-1/2},{j-1/2} flux limiter
				for (var j = this.nghost; j < this.nx+this.nghost-3; ++j) {
					for (var i = this.nghost; i < this.nx+this.nghost-3; ++i) {
						for (var side = 0; side < 2; ++side) {
							//dq = q_i,j - q_{{i,j}-dirs[side]}
							var dq = this.q[iq + 4 * (i + this.nx * j)]
								- this.q[iq + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))];
							if (Math.abs(dq) > 0) {
								if (this.ui[side + 2 * (i + (this.nx+1) * j)] >= 0) {
									this.r[iq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[iq + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))]
											- this.q[iq + 4 * (i - 2*dirs[side][0] + this.nx * (j - 2*dirs[side][1]))]) / dq;
								} else {
									this.r[iq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[iq + 4 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]))]
											- this.q[iq + 4 * (i + this.nx * j)]) / dq;
									
								}
							} else {
								this.r[iq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							}
						}
					}
				}
			
				//now for ghost boundaries
				for (var j = 0; j < this.nghost; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							this.r[iq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.r[iq + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.r[iq + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.r[iq + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
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
						/**/
						sideLengths[0] = dx;
						sideLengths[1] = dy;
					
						//flux calculation 
						for (var side = 0; side < 2; ++side) {
							var ir = iq + 4 * (side + 2 * (i + (this.nx+1) * j));
							//apply limiter
							var phi = this.fluxMethod(this.r[ir]);
							var iu = side + 2 * (i + (this.nx+1) * j);
							var iflux = iq + 4 * (side + 2 * (i + (this.nx+1) * j));
							var iqL = iq + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
							var iqR = iq + 4 * (i + this.nx * j);
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
							this.flux[iq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.flux[iq + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.flux[iq + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.flux[iq + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
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
							var ifluxR = iq + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
							var ifluxL = iq + 4 * (side + 2 * (i + (this.nx+1) * j));
							var df = this.flux[ifluxR] - this.flux[ifluxL];
							this.q[iq + 4 * (i + this.nx * j)] -= dt * df / dx;
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
						var iqL = ix - dirs[side][0] + this.nx * (iy - dirs[side][1]);
						var densityL = this.q[0 + iqL];
						var velocityXL = this.q[1 + iqL] / densityL;
						var velocityYL = this.q[2 + iqL] / densityL;
						var energyTotalL = this.q[3 + iqL] / densityL;
						var energyKinematicL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
						var energyThermalL = energyTotalL - energyKinematicL;
						var pressureL = (this.gamma - 1) * densityL * energyThermalL;
						var speedOfSoundL = Math.sqrt(this.gamma * prsesureL / densityL);
						var hTotalL = energyTotalL + pressureL / densityL;
						var roeWeightL = Math.sqrt(densityL);
						
						var iqR = ix + this.nx * iy;
						var densityR = this.q[0 + iqR];
						var velocityXR = this.q[1 + iqR] / densityR;
						var velocityYR = this.q[2 + iqR] / densityR;
						var energyTotalR = this.q[3 + iqR] / densityR;
						var energyKinematicR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
						var energyThermalR = energyTotalR - energyKinematicR;
						var pressureR = (this.gamma - 1) * densityR * energyThermalR;
						var speedOfSoundR = Math.sqrt(this.gamma * prsesureR / densityR);
						var hTotalR = energyTotalR + pressureR / densityR;
						var roeWeightR = Math.sqrt(densityR);
					
						var denom = roeWeightL + roeWeightR;
						var velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
						var velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
						var hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;

						buildEigenstate(
							 side + 2 * (ix + (this.nx+1) * iy),	//index into interface element.  from there you'll have to scale by cell size.  Thus manually recreating the automatic storage of C structures. JavaScript, why can't you be more like LuaJIT? 
							 this.interfaceMatrix,
							 this.interfaceEigenvalues,
							 this.interfaceEigenvectors,
							 this.interfaceEigenvectorsInverse,
							 velocityX, velocityY, hTotal, this.gamma,
							 dirs[side][0], dirs[side][1]);

						var rho = this.q[0 + 4 * (ix + this.nx * iy)];
						var u = this.q[1 + 4 * (ix + this.nx * iy)] / rho;
						var v = this.q[2 + 4 * (ix + this.nx * iy)] / rho; 
						var energyTotal = this.q[3 + 4 * (ix + this.nx * iy)] / rho; 
						var energyKinematic = .5 * (u * u + v * v);
						var energyThermal = energyTotal - energyKinematic;
						var speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
						var dx = this.xi[0 + 2 * (ix+1 + (this.nx+1) * iy)] - this.xi[0 + 2 * (ix + (this.nx+1) * iy)];
						var dy = this.xi[1 + 2 * (ix + (this.nx+1) * (iy+1))] - this.xi[1 + 2 * (ix + (this.nx+1) * iy)];
						//http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
						//var dum = 1 / (speedOfSound + Math.abs(u) / dx + Math.abs(v) / dy); 
						var dum = dx / (speedOfSound + Math.abs(u));
						if (mindum === undefined || dum < mindum) mindum = dum;
						var dum = dy / (speedOfSound + Math.abs(v));
						if (mindum === undefined || dum < mindum) mindum = dum;
					}
				}
			}
			return this.cfl * mindum;
	
		},
		advect : function(dt) {
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
		
		//f_{i-1/2},{j-1/2}: cell flux
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
					for (var jq = 0; jq < 4; ++jq) {
						for (var iq = 0; iq < 4; ++iq) {
							//initialize to identity matrix
							this.interfaceMatrix[iq + 4 * (jq + 4 * (side + 2 * (i + (this.nx+1) * j)))] = iq == jq ? 1 : 0;
							this.interfaceEigenvectors[iq + 4 * (jq + 4 * (side + 2 * (i + (this.nx+1) * j)))] = iq == jq ? 1 : 0;
							this.interfaceEigenvectorsInverse[iq + 4 * (jq + 4 * (side + 2 * (i + (this.nx+1) * j)))] = iq == jq ? 1 : 0;
						}
						this.interfaceEigenvalues[jq + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
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
		var ip = 0;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var rho = this.q[0 + iq];
				var u = this.q[1 + iq] / rho; 
				var v = this.q[2 + iq] / rho; 
				var energyTotal = this.q[3 + iq] / rho; 
				var energyKinematic = .5 * (u * u + v * v);
				var energyThermal = energyTotal - energyKinematic;
				this.pressure[ip] = (this.gamma - 1) * rho * energyThermal;
				++ip;
				iq += 4;
			}
		}
		
		//apply momentum diffusion = pressure
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var iq = 4 * (i + this.nx * j);
				for (var side = 0; side < 2; ++side) {
					var plusIndex = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					var minusIndex = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					this.q[1+side + iq] -= dt 
						* (this.pressure[plusIndex] - this.pressure[minusIndex]) 
						/ (this.x[side + 2 * plusIndex] - this.x[side + 2 * minusIndex]);
				}
			}
		}

		//apply work diffusion = momentum
		for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
				var iq = 4 * (i + this.nx * j);
				for (var side = 0; side < 2; ++side) {
					var iqR = 4 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]));
					var iqL = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
					//this is pulling the coordinate associated with the interface's direction
					//a more robust method would be to take both velocity components and dot them with the interface tangent
					var uR = this.q[1+side + iqR] / this.q[0 + iqR];
					var uL = this.q[1+side + iqL] / this.q[0 + iqL];
					var ipR = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					var ipL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					var ixR = side + 2 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]));
					var ixL = side + 2 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
					var dx = this.x[ixR] - this.x[ixL];
					this.q[3 + iq] -= dt * (this.pressure[ipR] * uR - this.pressure[ipL] * uL) / dx;
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
	$(window).resize(onresize);
	onresize();
	update();
});
