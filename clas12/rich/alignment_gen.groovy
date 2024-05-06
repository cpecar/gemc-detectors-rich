import org.jlab.geom.base.*;
import org.jlab.detector.geant4.v2.RICHGeant4Factory;
import org.jlab.detector.geom.RICH.RICHGeoFactory;
import org.jlab.detector.geom.RICH.RICHGeoCalibration;
import org.jlab.detector.geom.RICH.RICHGeoParameters;
import org.jlab.detector.geom.RICH.RICHLayer;
import org.jlab.detector.geom.RICH.RICHFrame;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.volume.*;
import org.jlab.geom.prim.*;
import java.lang.Math;

def variation = args[0];
int runnum = (args[1]).toInteger();

def outfile = "transformations_${variation}.txt"
ConstantsManager ccdb_geoFactory = new ConstantsManager();
ConstantsManager ccdb = new ConstantsManager();
String[] richTables = new String[]{
                        "/geometry/rich/setup",
                        "/geometry/rich/geo_parameter",
                        "/geometry/rich/module1/aerogel",
                        "/geometry/rich/module2/aerogel",
                        "/geometry/rich/module1/alignment",
                        "/geometry/rich/module2/alignment"
                     };
ccdb_geoFactory.setVariation("default");
ccdb_geoFactory.init(Arrays.asList(richTables));

ccdb.setVariation(variation);
ccdb.init(Arrays.asList(richTables));

System.out.println(ccdb);

// create RICH coatjava geometry objects and initialize
RICHGeoFactory geoFactory = new RICHGeoFactory(1, ccdb_geoFactory, runnum, true);
RICHGeoCalibration geoCal = new RICHGeoCalibration();
RICHGeoParameters geoPar = new RICHGeoParameters();

geoPar.load_CCDB(ccdb, runnum, 0, true);
geoCal.load_CCDB(ccdb, runnum, 0, geoPar);

// GEMC - lab frame
// RICH alignment ccdb - RICH frame, with rotations applied first
// after shifting by -1*barycenter of layer
// Get all values from ccdb, transform to correct lab frame angles/shifts, store in text file

// get layer frame definition

// align_Element goes like:
// if angle > 0:
//    translate shape in -1*frame barycenter (layerframe.bref())
//    3d angle into lab frame (from rich frame or layer frame)
//    rotate z, y, x
//    translate	shape in 1*frame barycenter (frame.bref())

// if shift > 0:
//    shift into lab frame (from rich frame?)
//    translate by shift_lab

// layer 201-204,401: all we need is global shift/rotation of plane
// just re-implementing method from RICHGeoFactory
def toLabFrame(Vector3D vec, RICHFrame frame){
    xref = frame.xref();
    yref = frame.yref();
    zref = frame.zref();
    return new Vector3D(         vec.z()*zref.x() + vec.y()*yref.x() + vec.x()*xref.x(),
                             vec.z()*zref.y() + vec.y()*yref.y() + vec.x()*xref.y(),
                             vec.z()*zref.z() + vec.y()*yref.z() + vec.x()*xref.z());
			     }

// invert rotation matrix given rotation around z, y, and x
def toInvertedR(alpha,beta,gamma){ // z, y, x
    def r11 = Math.cos(alpha) * Math.cos(beta)
    def r12 = Math.cos(alpha) * Math.sin(beta) * Math.sin(gamma) - Math.sin(alpha) * Math.cos(gamma)
    def r13 = Math.cos(alpha) * Math.sin(beta) * Math.cos(gamma) + Math.sin(alpha) * Math.sin(gamma)
    def r21 = Math.sin(alpha) * Math.cos(beta)
    def r22 = Math.sin(alpha) * Math.sin(beta) * Math.sin(gamma) + Math.cos(alpha) * Math.cos(gamma)
    def r23 = Math.sin(alpha) * Math.sin(beta) * Math.cos(gamma) - Math.cos(alpha) * Math.sin(gamma)
    def r31 = -Math.sin(beta)
    def r32 = Math.cos(beta) * Math.sin(gamma)
    def r33 = Math.cos(beta) * Math.cos(gamma)

    def rTranspose = [[r11, r21, r31], [r12, r22, r32], [r13, r23, r33]]

    def phi = Math.atan2(rTranspose[2][1], rTranspose[2][2]);
    def theta = -Math.asin(rTranspose[2][0]);
    def psi = Math.atan2(rTranspose[1][0], rTranspose[0][0]);
    
    return [phi,theta,psi];
}

sector = 4;

new File(outfile).withWriter('UTF-8') { writer ->

layers_global = [0,1,2,3,12,4,5,6,7,8,9,10] // layers where each component only gets global layer alignment
// AEROGEL+MAPMT LOOP: aligning each full aerogel/mapmt layer
for(int i = 0; i < layers_global.size(); i++){
	RICHLayer layer = geoFactory.get_Layer(4,layers_global[i]);
	
	System.out.println(layer.name());
	System.out.println("layer size: " + layer.size());
	def layerid = layer.id();
	
	RICHFrame lframe = layer.generate_LocalRef();
	RICHFrame rich_frame;
	if(sector == 4){
		rich_frame = geoFactory.richframes.get(0);
	}
	else{
		rich_frame = geoFactory.richframes.get(1);
	}
	Vector3D mothershift = geoCal.get_AlignShift(4,0,0);
	Vector3D motherangle = geoCal.get_AlignAngle(4,0,0);

	Vector3D lshift = geoCal.get_AlignShift(4,layerid+1,0);
	Vector3D langle = geoCal.get_AlignAngle(4,layerid+1,0);

	System.out.println("lshift: " + lshift);
	System.out.println("langle: " + langle);
	
	Vector3D lbary = lframe.bref();

	Vector3D langle_lab = toLabFrame(langle, lframe);
	Vector3D lshift_lab = toLabFrame(lshift, lframe);

	Vector3D motherangle_lab = toLabFrame(motherangle, rich_frame);
	Vector3D mothershift_lab = toLabFrame(mothershift, rich_frame);

	
	xang = langle_lab.x;
	yang = langle_lab.y;
	zang = langle_lab.z;
        println "l$layerid angles in lab frame: $xang $yang $zang";

	xshift = lshift_lab.x;
	yshift = lshift_lab.y;
	zshift = lshift_lab.z;
	println "l$layerid shift in lab frame: $xshift $yshift $zshift"

	baryx = lbary.x;
	baryy = lbary.y;
	baryz = lbary.z;	

	// RICHGeoFactory first shifts by vector of RICHbarycenter,
	// so full transformation is Rp + (Rv - v) (verify this)
	// R - rotation matrix, v - barycenter vector
	
	lbary_rotated = new Vector3D(lbary);
	lbary_rotated.rotateZ(zang);
	lbary_rotated.rotateY(yang);
	lbary_rotated.rotateX(xang);
	
	bary_shift_final = lbary.sub(lbary_rotated);
	xshift = bary_shift_final.x + lshift_lab.x;
	yshift = bary_shift_final.y + lshift_lab.y;
	zshift = bary_shift_final.z + lshift_lab.z;
	
	// G4PVPlacement expects inverse of rotation matrix
	newAngles = toInvertedR(zang,yang,xang);
	xang = newAngles[0]
	yang = newAngles[1]
	zang = newAngles[2]
	
	println "l$layerid angles after inversion: $xang $yang $zang";
	
	//String line = "$sector $layerid 0 $xang $yang $zang $xshift $yshift $zshift";
	String line = String.format("%d %d 0 %.8f %.8f %.8f %.8f %.8f %.8f ", sector, layerid, xang, yang, zang, xshift, yshift, zshift)

	writer.writeLine(line);
}

// spherical mirror alignment: individual component only
RICHLayer layer = geoFactory.get_Layer(4,11);
System.out.println(layer.name());
System.out.println("layer size: " + layer.size());
def layerid = layer.id();
RICHFrame lframe = layer.generate_LocalRef();

for(int i = 0; i < layer.size(); i++){
	RICHFrame cframe = layer.generate_LocalRef(i);
	Vector3D cshift = geoCal.get_AlignShift(sector,layerid+1,i+1);
        Vector3D cangle = geoCal.get_AlignAngle(sector,layerid+1,i+1);
	
	// just cframe instead of lframe
 	Vector3D cbary = cframe.bref();

        Vector3D cangle_lab = toLabFrame(cangle, cframe);
        Vector3D cshift_lab = toLabFrame(cshift, cframe);

        xang = cangle_lab.x;
        yang = cangle_lab.y;
        zang = cangle_lab.z;
        println "l$layerid angles in lab frame: $xang $yang $zang"
        xshift = cshift_lab.x;
        yshift = cshift_lab.y;
        zshift = cshift_lab.z;
        println "l$layerid shift in lab frame: $xshift $yshift $zshift"

        baryx = cbary.x;
        baryy = cbary.y;
        baryz = cbary.z;

        // RICHGeoFactory first shifts by vector of RICHbarycenter,
        // so full transformation is Rp + (Rv - v) (verify this)
        // R - rotation matrix, v - barycenter vector

        cbary_rotated = new Vector3D(cbary);
        cbary_rotated.rotateZ(zang);
        cbary_rotated.rotateY(yang);
        cbary_rotated.rotateX(xang);

        bary_shift_final = cbary.sub(cbary_rotated);
        xshift = bary_shift_final.x + cshift_lab.x;
        yshift = bary_shift_final.y + cshift_lab.y;
        zshift = bary_shift_final.z + cshift_lab.z;		
	
	// G4PVPlacement expects inverse of rotation matrix
	newAngles = toInvertedR(zang,yang,xang);
	xang = newAngles[0]
	yang = newAngles[1]
	zang = newAngles[2]
	
        String line = "$sector $layerid $i $xang $yang $zang $xshift $yshift $zshift";
        writer.writeLine(line);


}

}