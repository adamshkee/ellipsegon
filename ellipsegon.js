// returns the array of points for a leaflet polygon that correspond to an ellipse
// with the given parameters:
// lat_deg  - lattitude of centre point in degrees
// lon_deg  - longitude of centre point in degrees
// maj_m    - semi-major-axis on metres
// min_m    - semi-minor-axis in metres
// tilt_deg - angle counter-clockwise from west to the semi-major-axis in degrees
// npoints  - number of points to use in the ellipse
function getEllipsePoints(lat_deg, lon_deg, maj_m, min_m, tilt_deg, npoints) {
	let points = [];

	for (let i = 0; i < npoints; i++) {
		let angle_deg = i * 360.0 / npoints;
		let radius_m = getEllipseRadiusAtAngle(maj_m, min_m, tilt_deg, angle_deg);
		let latLon = getEndPoint(lat_deg, lon_deg, angle_deg, radius_m);
		points.push(latLon);
	}

	return points;
}

// Given an ellipse with the given parameters, return radius in at a given angle in degrees
// semimajor and semi minor axis (maj, min) are assumed to be in the same units, and the
// return value will be these units as well.
function getEllipseRadiusAtAngle(maj, min, tilt_deg, angle_deg) {
	let ab = maj * min;
	let a2 = Math.pow(maj, 2);
	let b2 = Math.pow(min, 2);

	let thetaPrime_rad = (90 - angle_deg - tilt_deg) * Math.PI / 180.0;
	let sinThetaPrime = Math.sin(thetaPrime_rad);
	let cosThetaPrime = Math.cos(thetaPrime_rad);
	let sin2ThetaPrime = Math.pow(sinThetaPrime, 2);
	let cos2ThetaPrime = Math.pow(cosThetaPrime, 2);

	let r = ab / Math.sqrt(a2 * sin2ThetaPrime + b2 * cos2ThetaPrime);

	return r;
}

// Given starting lat/lon in degrees, initial bearing in degrees and distance in metres
// returns end lat/lon in degrees as an array [lat, lon]
function getEndPoint(lat_deg, lon_deg, bearing_deg, distance_m) {
// courtesy of https://www.movable-type.co.uk/scripts/latlong-vincenty.html
	let a = 6378137; //Equitorial Radius metres
	let b = 6356752.3; //Polar Radius metres
	let f = (a-b)/a;
	
	let lat_rad = lat_deg * Math.PI / 180.0;
	let lon_rad = lon_deg * Math.PI / 180.0;	

	let bearing_rad = bearing_deg * Math.PI / 180.0;

	let sin_bearing = Math.sin(bearing_rad);
	let cos_bearing = Math.cos(bearing_rad);

	let tanU1 = (1-f) * Math.tan(lat_rad);
	let cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1));
	let sinU1 = tanU1 * cosU1;

	let sigma1 = Math.atan2(tanU1, cos_bearing);
	let sina = cosU1 * sin_bearing;
	let cosSqa = 1 - sina*sina;
	let uSq = cosSqa * (a*a - b*b) / (b*b);
	let A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
	let B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

	let sigma = distance_m / (b*A), sigma_prime;
	do {
	    var cos2sigmaM = Math.cos(2*sigma1 + sigma);
	    var sinsigma = Math.sin(sigma);
	    var cossigma = Math.cos(sigma);
	    let deltasigma = B*sinsigma*(cos2sigmaM+B/4*(cossigma*(-1+2*cos2sigmaM*cos2sigmaM)-
		B/6*cos2sigmaM*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmaM*cos2sigmaM)));
	    sigma_prime = sigma;
	    sigma = distance_m / (b*A) + deltasigma;
	} while (Math.abs(sigma-sigma_prime) > 1e-12);

	let tmp = sinU1*sinsigma - cosU1*cossigma*cos_bearing;
	let lat2_rad = Math.atan2(sinU1*cossigma + cosU1*sinsigma*cos_bearing, (1-f)*Math.sqrt(sina*sina + tmp*tmp));
	let lambda = Math.atan2(sinsigma*sin_bearing, cosU1*cossigma - sinU1*sinsigma*cos_bearing);
	let C = f/16*cosSqa*(4+f*(4-3*cosSqa));
	let L = lambda - (1-C) * f * sina * (sigma + C*sinsigma*(cos2sigmaM+C*cossigma*(-1+2*cos2sigmaM*cos2sigmaM)));
	let lon2_rad = (lon_rad+L+3*Math.PI)%(2*Math.PI) - Math.PI;  // normalise to -180...+180

	let revAz = Math.atan2(sina, -tmp);
	
	return [lat2_rad * 180.0 / Math.PI, lon2_rad * 180.0 / Math.PI];
}