const fluxMethods = {
	'donor cell' : r => { return 0; },
	'Lax-Wendroff' : r => { return 1; },
	
	//these two are no good with the Godunov (Riemann) solver
	'Beam-Warming' : r => { return r; },
	'Fromm' : r => { return .5 * (1 + r); },

	//Wikipedia
	CHARM : r => { if (r < 0) return 0; return r*(3*r+1)/((r+1)*(r+1)); },
	HCUS : r => { return Math.max(0, 1.5 * (r + Math.abs(r)) / (r + 2) ); },
	HQUICK : r => { return Math.max(0, 2 * (r + Math.abs(r)) / (r + 3) ); },
	Koren : r => { return Math.max(0, Math.min(2*r, (1 + 2*r)/3 ,2) ); },
	minmod : r => { return Math.max(0, Math.min(r,1) ); },
	Oshker : r => { return Math.max(0, Math.min(r,1.5) ); },	//replace 1.5 with 1 <= beta <= 2	
	ospre : r => { return .5 * (r*r + r) / (r*r + r + 1); },
	smart : r => { return Math.max(0, Math.min(2 * r, .25 + .75 * r, 4)); },
	Sweby : r => { return Math.max(0, Math.min(1.5 * r, 1), Math.min(r, 1.5)); },	//replace 1.5 with 1 <= beta <= 2
	UMIST : r => { return Math.max(0, Math.min(2*r, .75 + .25*r, .25 + .75*r, 2)); },	
	'van Albada 1' : r => { return (r * r + r) / (r * r + 1); },
	'van Albada 2' : r => { return 2 * r / (r * r + 1); },
	
	'van Leer' : r => { return (r + Math.abs(r)) / (1 + Math.abs(r)); },
	'monotonized central' : r => { return Math.max(0, Math.min(2, .5 * (1 + r), 2 * r)); },
	superbee : r => { return Math.max(0,Math.min(1,2*r),Math.min(2,r)); },

	'Barth-Jespersen' : r => { return .5 * (r + 1) * Math.min(1, 4*r/(r+1), 4/(r+1)); }
};
export {fluxMethods};
