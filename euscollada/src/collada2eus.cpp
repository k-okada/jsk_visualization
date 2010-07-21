#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#include "dae.h"
#include "dom/domCOLLADA.h"

#include <fstream>
#include "yaml-cpp/yaml.h"

#include <boost/foreach.hpp>

daeDocument *g_document;
DAE* g_dae = NULL;

vector<pair<string, string> > all_link_names;

// returns max offsset value
unsigned int getMaxOffset( domInput_local_offset_Array &input_array )
{
  unsigned int maxOffset = 0;
  for ( unsigned int i = 0; i < input_array.getCount(); i++ ) {
    if ( input_array[i]->getOffset() > maxOffset ) {
      maxOffset = (unsigned int)input_array[i]->getOffset();
    }
  }
  return maxOffset;
}

void writeTriangle(FILE *fp, domGeometry *thisGeometry) {

  // get mesh
  domMesh *thisMesh = thisGeometry->getMesh();
  int triangleElementCount = (int)(thisMesh->getTriangles_array().getCount());

  fprintf(fp, "(defclass %s\n", thisGeometry->getId());
  fprintf(fp, "  :super bodyset-link\n");
  fprintf(fp, "  :slots ())\n");
  fprintf(fp, "(defmethod %s\n", thisGeometry->getId());
  fprintf(fp, "  (:init (&key (name))\n");
  fprintf(fp, "         (send-super :init (make-cascoords) :bodies (list (make-cube 10 10 10)) :name name))\n");
  fprintf(fp, "  (:draw (vwr)\n");
  fprintf(fp, "   (let ((mat (send (send self :worldcoords) :4x4))\n");
  fprintf(fp, "#+:jsk\n");
  fprintf(fp, "         (glcon (cdr (assq (sys:thread-self)\n");
  fprintf(fp, "                           ((send vwr :viewsurface) . gl::glcon))))\n");
  fprintf(fp, "#-:jsk\n");
  fprintf(fp, "         (glcon ((send vwr :viewsurface) . gl::glcon))\n");
  fprintf(fp, "         newlis)\n");
  fprintf(fp, "#+:jsk\n");
  fprintf(fp, "     (sys::mutex-lock gl::*opengl-lock*)\n");
  fprintf(fp, "     (send vwr :viewsurface :makecurrent)\n");
  fprintf(fp, "     (gl::glPushAttrib gl::GL_ALL_ATTRIB_BITS)\n");
  fprintf(fp, "     (gl::glPushMatrix)\n");
  fprintf(fp, "     (gl::glMultMatrixf (array-entity (transpose mat gl::*temp-matrix*)))\n");
  fprintf(fp, "     (if (setq newlis (cdr (assq glcon (get self :GL-DISPLAYLIST-ID))))\n");
  fprintf(fp, "         (progn\n");
  fprintf(fp, "           (gl::glCallList newlis)\n");
  fprintf(fp, "           )\n");
  fprintf(fp, "       (progn\n");
  fprintf(fp, "         (setq newlis (gl::glGenLists 1))\n");
  fprintf(fp, "         (gl::glNewList newlis gl::GL_COMPILE)\n");
  // Triangulation
  // based on http://www.shader.jp/xoops/html/modules/mydownloads/singlefile.php?cid=5&lid=6
  for(int currentTriangle=0;currentTriangle<triangleElementCount;currentTriangle++)
    {
      domTriangles* thisTriangles = thisMesh->getTriangles_array().get(0);
      domMaterial* thisMaterial = daeSafeCast<domMaterial>(g_dae->getDatabase()->idLookup(string(thisGeometry->getId())+string(".mat"), g_document));
      domInstance_effect* thisInstanceEffect = thisMaterial->getInstance_effect();
      domEffect* thisEffect = daeSafeCast<domEffect>(g_dae->getDatabase()->idLookup(thisInstanceEffect->getUrl().id(),g_document));

      fprintf(fp, "         (gl::glColor3fv (float-vector 0.1 0.1 0.1))\n");
      fprintf(fp, "         (gl::glMaterialfv gl::GL_FRONT gl::GL_AMBIENT (float-vector %f %f %f %f))\n",
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getAmbient()->getColor()->getValue()[0],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getAmbient()->getColor()->getValue()[1],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getAmbient()->getColor()->getValue()[2],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getAmbient()->getColor()->getValue()[3]);
      fprintf(fp, "         (gl::glMaterialfv gl::GL_FRONT gl::GL_DIFFUSE (float-vector %f %f %f %f))\n",
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getDiffuse()->getColor()->getValue()[0],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getDiffuse()->getColor()->getValue()[1],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getDiffuse()->getColor()->getValue()[2],
              thisEffect->getFx_profile_array()[0]->getProfile_COMMON()->getTechnique()->getPhong()->getDiffuse()->getColor()->getValue()[3]);


      int numberOfInputs = (int)getMaxOffset(thisTriangles->getInput_array()) +1;// offset
      int numberOfTriangles = (int)(thisTriangles->getP()->getValue().getCount() / numberOfInputs);	// elements

      // offset of index
      unsigned int offset = 0;
      int texoffset = -255, noroffset = -255;
      for(unsigned int i=0;i<thisTriangles->getInput_array().getCount();i++)
        {
          if(strcmp(thisTriangles->getInput_array()[i]->getSemantic(), "VERTEX")==0)
            offset = thisTriangles->getInput_array()[i]->getOffset();
          if(strcmp(thisTriangles->getInput_array()[i]->getSemantic(), "TEXCOORD")==0)
            texoffset = thisTriangles->getInput_array()[i]->getOffset();
          if(strcmp(thisTriangles->getInput_array()[i]->getSemantic(), "NORMAL")==0)
            noroffset = thisTriangles->getInput_array()[i]->getOffset();
        }
      fprintf(fp, "         (gl::glBegin gl::GL_TRIANGLES)\n");
      for(int i=0;i<numberOfTriangles;i++)
        {
          int index = thisTriangles->getP()->getValue().get(i*numberOfInputs+offset);
          int sourceElements = thisMesh->getSource_array().getCount();

          // normal
          if ( sourceElements  > 1 ) {
            if(noroffset==-255)
              {
                // normal vectur shares same index of vertices
                fprintf(fp, "         (gl::glNormal3fv (float-vector %f %f %f))\n",
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(index*3),
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(index*3+1),
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(index*3+2)
                        );
              }
            else
              {
                // index normal vector is indicated in <p></p>
                int norindex = thisTriangles->getP()->getValue().get(i*numberOfInputs+noroffset);
                fprintf(fp, "         (gl::glNormal3f (float-vector %f %f %f))\n",
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(norindex*3),
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(norindex*3+1),
                        thisMesh->getSource_array()[1]->getFloat_array()->getValue().get(norindex*3+2)
                        );
              }
          } else {
            if ( i % 3 == 0 ) {
              int normal_index_0 = thisTriangles->getP()->getValue().get((i+0)*numberOfInputs+offset);
              int normal_index_1 = thisTriangles->getP()->getValue().get((i+1)*numberOfInputs+offset);
              int normal_index_2 = thisTriangles->getP()->getValue().get((i+2)*numberOfInputs+offset);
              float a0 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_2)*3+0) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+0);
              float a1 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_2)*3+1) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+1);
              float a2 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_2)*3+2) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+2);
              float b0 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_1)*3+0) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+0);
              float b1 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_1)*3+1) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+1);
              float b2 =
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_1)*3+2) -
                thisMesh->getSource_array()[0]->getFloat_array()->getValue().get((normal_index_0)*3+2);
              float aa = sqrt(a0*a0+a1*a1+a2*a2); a0 /= aa; a1 /= aa; a2 /= aa;
              float bb = sqrt(b0*b0+b1*b1+b2*b2); b0 /= bb; b1 /= bb; b2 /= bb;

              fprintf(fp, "         (gl::glNormal3fv (float-vector %f %f %f))\n",
                      (a1*b2 - a2*b1), (a2*b0 - a0*b2), (a0*b1 - a1*b0));
            }
          }

          if ( sourceElements  > 2 ) {
            // texture coordinates
            if(texoffset!=-255)
              {
                int texindex = thisTriangles->getP()->getValue().get(i*numberOfInputs+texoffset);
                fprintf(fp, "         (gl::glTexCoords2f %f %f)\n",
                        thisMesh->getSource_array()[2]->getFloat_array()->getValue().get(texindex*2),
                        thisMesh->getSource_array()[2]->getFloat_array()->getValue().get(texindex*2+1)
                        );
              }
          }

          // vertex vector
          fprintf(fp, "         (gl::glVertex3fv (float-vector %f %f %f))\n",
                  1000*thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3),
                  1000*thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3+1),
                  1000*thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3+2)
                  );
        }
      fprintf(fp, "         (gl::glEnd)\n");
      fprintf(fp, "         (gl::glEndList)\n");
      fprintf(fp, "         (setf (get self :GL-DISPLAYLIST-ID)\n");
      fprintf(fp, "               (cons (cons glcon newlis)\n");
      fprintf(fp, "                     (get self :GL-DISPLAYLIST-ID)))\n");
      fprintf(fp, "         (setq newlis nil)))\n");
      fprintf(fp, "     (gl::glPopMatrix)\n");
      fprintf(fp, "     (gl::glPopAttrib)\n");
      fprintf(fp, "#+:jsk\n");
      fprintf(fp, "     (sys::mutex-unlock gl::*opengl-lock*)\n");
      fprintf(fp, "     (unless newlis (send self :draw vwr))\n");
      fprintf(fp, "     )))\n");
      fprintf(fp, "  )\n");
    }

  // Polylist
  int polylistElementCount = (int)(thisMesh->getPolylist_array().getCount());
  //Polygons
  int polygonesElementCount = (int)(thisMesh->getPolygons_array().getCount());
  assert(polylistElementCount==0);
  assert(polygonesElementCount==0);
}

void writeGeometry(FILE *fp, daeDatabase *thisDatabase) {
  // number of geometry
  int geometryElementCount =  thisDatabase->getElementCount(NULL, "geometry", NULL);
  for(int currentGeometry=0;currentGeometry<geometryElementCount;currentGeometry++) {
    // get current geometry
    domGeometry *thisGeometry;
    thisDatabase->getElement((daeElement**)&thisGeometry, currentGeometry, NULL, "geometry");

    fprintf(stderr, "geometry %d id:%s\n",
            currentGeometry, thisGeometry->getId());

    // write geometry information
    writeTriangle(fp, thisGeometry);
  }
}

//
const char *findGeometryFromLinkName(const char *linkName) {
  int nodeElementCount;
  nodeElementCount = g_dae->getDatabase()->getElementCount(NULL, "node", NULL);
  for(int currentNode=0;currentNode<nodeElementCount;currentNode++) {
    domNode *thisNode;
    g_dae->getDatabase()->getElement((daeElement**)&thisNode, currentNode, NULL, "node");
    if (strcmp(linkName, thisNode->getName())==0) {
      while ( thisNode->getInstance_geometry_array().getCount()==0 &&
              thisNode->getNode_array().getCount() != 0 ) {
        thisNode = thisNode->getNode_array()[0];
      }
      int geometryElementCount = thisNode->getInstance_geometry_array().getCount();
      for(int currentGeometry=0;currentGeometry<geometryElementCount;currentGeometry++) {
        domInstance_geometry *thisGeometry = thisNode->getInstance_geometry_array()[currentGeometry];
        return thisGeometry->getUrl().id().c_str();
      }
    }
  }
  return NULL;
}

domJoint *findJointFromName(const char *jointName) {
  int jointElementCount;
  jointElementCount = g_dae->getDatabase()->getElementCount(NULL, "joint", NULL);
  for(int currentJoint=0;currentJoint<jointElementCount;currentJoint++) {
    domJoint *thisJoint;
    // get current geometry
    g_dae->getDatabase()->getElement((daeElement**)&thisJoint, currentJoint, NULL, "joint");
    string jointSid_str = string(((domKinematics_model *)(thisJoint->getParentElement()->getParentElement()))->getId())+string("/")+string(thisJoint->getSid());
    if ( jointSid_str.compare(jointName) == 0 ||
         strcmp(thisJoint->getName(),jointName) == 0 ) {
      return thisJoint;
    }
  }
  fprintf(stderr, "Counld not found joint %s\n", jointName);
  exit(1);
}

domLink *findNonFixChildLink(domLink *thisLink) {
  int attachmentCount = (int)(thisLink->getAttachment_full_array().getCount());
  for(int currentAttachment=0;currentAttachment<attachmentCount;currentAttachment++){
    domJoint *thisJoint = findJointFromName(thisLink->getAttachment_full_array()[currentAttachment]->getJoint());
    domAxis_constraint_Array jointAxis_array;
    int jointCount = 0;
    if ( thisJoint->getPrismatic_array().getCount() > 0 ) {
      jointAxis_array = thisJoint->getPrismatic_array();
      jointCount = thisJoint->getPrismatic_array().getCount();
    } else if ( thisJoint->getRevolute_array().getCount() > 0 ) {
      jointAxis_array = thisJoint->getRevolute_array();
      jointCount = thisJoint->getRevolute_array().getCount();
    }
    //
    if (jointAxis_array[0]->getLimits() &&
        jointAxis_array[0]->getLimits()->getMin()->getValue() == 0 &&
        jointAxis_array[0]->getLimits()->getMax()->getValue() == 0 ) {
      domLink *tmp_thisLink = findNonFixChildLink(thisLink->getAttachment_full_array()[currentAttachment]->getLink());
      if ( tmp_thisLink ) return tmp_thisLink;
    } else {
      return thisLink;
    }
  }
  return NULL;
}

domLink *findChildLinkFromJointName(const char *jointName) {
  domJoint* thisJoint = findJointFromName(jointName);
  string jointSid_str = string(((domKinematics_model *)(thisJoint->getParentElement()->getParentElement()))->getId())+string("/")+string(thisJoint->getSid());
  int linkElementCount = g_dae->getDatabase()->getElementCount(NULL, "link", NULL);
  for(int currentLink=0;currentLink<linkElementCount;currentLink++) {
    domLink *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    for(int currentAttachment=0;currentAttachment<(int)(thisLink->getAttachment_full_array().getCount());currentAttachment++){
      if ( jointSid_str.compare(thisLink->getAttachment_full_array()[currentAttachment]->getJoint()) == 0 ) {
        domLink *childLink = thisLink->getAttachment_full_array()[currentAttachment]->getLink();
        domLink *tmp_childLink = findNonFixChildLink(childLink);
        if (tmp_childLink) childLink = tmp_childLink;
        return childLink;
      }
    }
  }
  fprintf(stderr, "Counld not found joint %s\n", jointName);
  exit(1);
}

const char* findLinkName(const char *link_name) {
  for(vector<pair<string, string> >::iterator it=all_link_names.begin();it!=all_link_names.end();it++){
    if ( it->first.compare(link_name) == 0 ) {
      return it->second.c_str();
    }
  }
  return link_name;
}

// write euslisp joint instance from jointSid, parentLink, childLink
void writeJoint(FILE *fp, const char *jointSid, domLink *parentLink, domLink *childLink) {
  //
  domLink *tmp_childLink = findNonFixChildLink(childLink);
  if (tmp_childLink) childLink = tmp_childLink;

  // get number of joints
  int jointElementCount;
  jointElementCount = g_dae->getDatabase()->getElementCount(NULL, "joint", NULL);
  domJoint *thisJoint = findJointFromName(jointSid);
  fprintf(fp, "     (setq %s\n", thisJoint->getName());
  fprintf(fp, "           (instance %s :init\n",
          (thisJoint->getPrismatic_array().getCount()>0)?"linear-joint":"rotational-joint");
  fprintf(fp, "                     :name %s\n", findLinkName(thisJoint->getName()));
  fprintf(fp, "                     :parent-link %s :child-link %s\n", parentLink->getName(), childLink->getName());
  domAxis_constraint_Array jointAxis_array;
  int jointCount;
  float axis[3], min = -360, max = 360;
  if ( thisJoint->getPrismatic_array().getCount() > 0 ) {
    jointAxis_array = thisJoint->getPrismatic_array();
    if ( jointAxis_array[0]->getLimits() ) {
      min = 1000*jointAxis_array[0]->getLimits()->getMin()->getValue();
      max = 1000*jointAxis_array[0]->getLimits()->getMax()->getValue();
    }
  } else if ( thisJoint->getRevolute_array().getCount() > 0 ) {
    jointAxis_array = thisJoint->getRevolute_array();
    jointCount = thisJoint->getRevolute_array().getCount();
    if ( jointAxis_array[0]->getLimits() ) {
      min = jointAxis_array[0]->getLimits()->getMin()->getValue();
      max = jointAxis_array[0]->getLimits()->getMax()->getValue();
    }
  }
  axis[0] = jointAxis_array[0]->getAxis()->getValue()[0];
  axis[1] = jointAxis_array[0]->getAxis()->getValue()[1];
  axis[2] = jointAxis_array[0]->getAxis()->getValue()[2];
  fprintf(fp, "                     :axis (float-vector %f %f %f)\n", axis[0], axis[1], axis[2]);
  fprintf(fp, "                     :min %f :max %f\n", min, max);
  fprintf(fp, "                     ))\n");
}

void writeLink(FILE *fp, domLink::domAttachment_full_Array thisAttachmentArray) {

  int attachmentArrayCount = thisAttachmentArray.getCount();
  for(int currentAttachment=0;currentAttachment<attachmentArrayCount;currentAttachment++){
    domLinkRef thisLink = thisAttachmentArray[currentAttachment]->getLink();
    // link

    writeLink(fp, thisLink->getAttachment_full_array());

    // link
    fprintf(stderr, "link id:%s name:%s attachment:%d\n",
            thisLink->getSid(), thisLink->getName(), attachmentArrayCount);

    fprintf(fp, "\n");
    fprintf(fp, "     ;; define bodyset-link for %s\n", thisLink->getName());
    const char *geometry_name;
    if ( (geometry_name = findGeometryFromLinkName(thisLink->getName())) != NULL ) {
      fprintf(fp, "     (setq %s (instance %s :init :name \"%s\"))\n", thisLink->getName(), geometry_name, thisLink->getName());
    } else {
      fprintf(fp, "     (setq %s (instance bodyset-link :init (make-cascoords) :bodies (list (make-cube 1 1 1)) :name \"%s\"))\n",
              thisLink->getName(), thisLink->getName());
    }

    // assoc
    for(int currentAttachment2=0;currentAttachment2<(int)(thisLink->getAttachment_full_array().getCount());currentAttachment2++){
      fprintf(fp, "     (send %s :assoc %s)\n",
              thisLink->getName(),
              thisLink->getAttachment_full_array()[currentAttachment2]->getLink()->getName());
      writeJoint(fp, thisLink->getAttachment_full_array()[currentAttachment2]->getJoint(),
                 thisLink,
                 thisLink->getAttachment_full_array()[currentAttachment2]->getLink()
                 );
    }

    // translate
    int translateCount = thisAttachmentArray[currentAttachment]->getTranslate_array().getCount();
    for(int currentTranslate=0;currentTranslate<translateCount;currentTranslate++){
      domTranslateRef thisTranslate = thisAttachmentArray[currentAttachment]->getTranslate_array()[currentTranslate];
      fprintf(fp, "     (send %s :translate (float-vector %f %f %f))\n",
              thisLink->getName(),
              1000*thisTranslate->getValue()[0],
              1000*thisTranslate->getValue()[1],
              1000*thisTranslate->getValue()[2]);
    }
    // rotate
    int rotateCount = thisAttachmentArray[currentAttachment]->getRotate_array().getCount();
    for(int currentRotate=0;currentRotate<rotateCount;currentRotate++){
      domRotateRef thisRotate = thisAttachmentArray[currentAttachment]->getRotate_array()[currentRotate];
      fprintf(fp, "     (send %s :rotate %f (float-vector %f %f %f))\n",
              thisLink->getName(),
              (thisRotate->getValue()[3])*M_PI/180.0,
              thisRotate->getValue()[0],
              thisRotate->getValue()[1],
              thisRotate->getValue()[2]);
    }
  }
}

void writeKinematics(FILE *fp, domLinkRef thisLink) {
  int attachmentArrayCount = thisLink->getAttachment_full_array().getCount();
  writeLink(fp, thisLink->getAttachment_full_array());

  fprintf(stderr, "link id:%s name:%s attachment:%d\n",
          thisLink->getSid(), thisLink->getName(), attachmentArrayCount);
  fprintf(fp, "     ;; define bodyset-link for %s\n", thisLink->getName());
  fprintf(fp, "     (setq %s (instance bodyset-link :init (make-cascoords) :bodies (list (make-cube 10 10 10))))\n",
          thisLink->getName());
  for(int currentAttachment2=0;currentAttachment2<(int)(thisLink->getAttachment_full_array().getCount());currentAttachment2++){
    fprintf(fp, "     (send %s :assoc %s)\n",
            thisLink->getName(),
            thisLink->getAttachment_full_array()[currentAttachment2]->getLink()->getName());
    writeJoint(fp, thisLink->getAttachment_full_array()[currentAttachment2]->getJoint(),
               thisLink,
               thisLink->getAttachment_full_array()[currentAttachment2]->getLink()
               );
  }
  fprintf(fp, "     (send self :assoc %s)\n", thisLink->getName());
  fprintf(fp, "\n");
}

int main(int argc, char* argv[]){
  FILE *output_fp;
  char *input_filename, *yaml_filename, *output_filename;
  if (argc!=4) {
    fprintf(stderr, "Usage: %s <input dae filename> <input yaml filename> <output dae filename>\n",argv[0]);
    exit(-1);
  }

  input_filename  = argv[1];
  yaml_filename   = argv[2];
  output_filename = argv[3];
  output_fp = fopen(output_filename,"w");
  if ( output_fp == NULL ) {
    exit(-1);
  }
  fprintf(stderr, "Convert %s to %s\n", input_filename, output_filename);

  // init COLLADA
  g_dae = new DAE();
  int iRet = g_dae->load(input_filename);
  if ( iRet != DAE_OK ) {
    exit(1);
  }

  if ( g_dae->getDocCount() != 1 ) {
    fprintf(stderr, "Number of documnet is not 1\n");
    exit(1);
  }
  g_document = g_dae->getDoc(0);

  // read yaml
  vector<string> larm_joint_names, rarm_joint_names, lleg_joint_names, rleg_joint_names, head_joint_names, torso_joint_names;
  vector<string> larm_link_names, rarm_link_names, lleg_link_names, rleg_link_names, head_link_names, torso_link_names;

  typedef pair<vector<string> &, vector<string> & > link_joint;
  typedef pair<string, link_joint > link_joint_pair;
  link_joint_pair limbs[]
    = {  link_joint_pair("larm", link_joint(larm_link_names, larm_joint_names)),
         link_joint_pair("rarm", link_joint(rarm_link_names, rarm_joint_names)),
         link_joint_pair("lleg", link_joint(lleg_link_names, lleg_joint_names)),
         link_joint_pair("rleg", link_joint(rleg_link_names, rleg_joint_names)),
         link_joint_pair("head", link_joint(head_link_names, head_joint_names)),
         link_joint_pair("torso", link_joint(torso_link_names, torso_joint_names)) };

  ifstream fin(yaml_filename);
  if (fin.fail()) {
    fprintf(stderr, "Could not open %s.", yaml_filename);
    exit(-1);
  }
  YAML::Parser parser(fin);
  YAML::Node doc;
  parser.GetNextDocument(doc);

  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string>& link_names = limb.second.first;
    vector<string>& joint_names = limb.second.second;
    try {
      const YAML::Node& limb_doc = doc[limb_name];
      for(unsigned int i = 0; i < limb_doc.size(); i++) {
        const YAML::Node& n = limb_doc[i];
        for(YAML::Iterator it=n.begin();it!=n.end();it++) {
          string key, value; it.first() >> key; it.second() >> value;
          cerr << key << "/" << value << endl;
          joint_names.push_back(key);
          link_names.push_back(findChildLinkFromJointName(key.c_str())->getName());
          all_link_names.push_back(pair<string, string>(key, value));
        }
      }
    } catch(YAML::RepresentationException& e) {
    }
  }

  // get number of kinmatics
  int visualSceneCount;
  visualSceneCount = g_dae->getDatabase()->getElementCount(NULL, "visual_scene", NULL);
  fprintf(stderr, "Number of Visual Scene %d (= 1)\n", visualSceneCount); // this shoule be 1
  domVisual_scene *thisVisualscene;
  g_dae->getDatabase()->getElement((daeElement**)&thisKinematics, 0, NULL, "visual_scene");
  int nodeCount = thisVisualscene->getNode_array().getCount();
  fprintf(stderr, "Number of Nodes %d (= 1)\n", nodeCount); // this shoule be 1
  domNode* = thisVisualscene->getNode_array()[0];

  fprintf(stderr, "Visual_scene %s\n", thisNode-getName());
  exit(1);

  fprintf(output_fp, ";; DO NOT EDIT THIS FILE\n");
  fprintf(output_fp, ";;\n");
  fprintf(output_fp, ";; this file is automatically generated by\n");
  fprintf(output_fp, ";; %s $ ", get_current_dir_name());for(int i=0;i<argc;i++) fprintf(output_fp, "%s ", argv[i]); fprintf(output_fp, "\n");
  fprintf(output_fp, ";;\n");
  fprintf(output_fp, "\n");
  fprintf(output_fp, "\n");
  fprintf(output_fp, "\n");
#if 0
  fprintf(output_fp, "(defclass %s-robot\n", thisKinematics->getName());
  fprintf(output_fp, "  :super robot-model\n");
  fprintf(output_fp, "  :slots ())\n");
  fprintf(output_fp, "(defmethod %s-robot\n", thisKinematics->getName());
  fprintf(output_fp, "  (:init\n");
  fprintf(output_fp, "   (&rest args)\n");
  fprintf(output_fp, "   (let (");

  // all joint and link name
  for(int currentJoint=0;currentJoint<(int)(g_dae->getDatabase()->getElementCount(NULL, "joint", NULL));currentJoint++) {
    domJoint *thisJoint;
    g_dae->getDatabase()->getElement((daeElement**)&thisJoint, currentJoint, NULL, "joint");
    fprintf(output_fp, "%s ", thisJoint->getName());
  }
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "link", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    fprintf(output_fp, "%s ", thisLink->getName());
  }
  fprintf(output_fp, ")\n");

  // send super :init
  fprintf(output_fp, "     (send-super* :init args)\n");
  fprintf(output_fp, "\n");

  // write kinemtaics
  writeKinematics(output_fp, thisKinematics->getTechnique_common()->getLink_array()[0]);
#endif
  // end-coords

  fprintf(output_fp, "     ;; end coords\n");
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string> link_names = limb.second.first;

    if (link_names.size()>0) {
      fprintf(output_fp, "     (setq %s-end-coords (make-cascoords :coords (send %s :copy-worldcoords)))\n", limb_name.c_str(), link_names.back().c_str());
      try {
        const YAML::Node& n = doc[limb_name+"-end-coords"]["translate"];
        double value;
        fprintf(output_fp, "     (send %s-end-coords :translate (float-vector", limb_name.c_str());
        for(unsigned int i=0;i<3;i++) { n[i]>>value; fprintf(output_fp, " %f", 1000*value);}
        fprintf(output_fp, "))\n");
      } catch(YAML::RepresentationException& e) {
      }
      try {
        const YAML::Node& n = doc[limb_name+"-end-coords"]["rotate"];
        double value;
        fprintf(output_fp, "     (send %s-end-coords :rotate", limb_name.c_str());
        for(unsigned int i=3;i<4;i++) { n[i]>>value; fprintf(output_fp, " %f", M_PI/180*value);}
        fprintf(output_fp, " (float-vector");
        for(unsigned int i=0;i<3;i++) { n[i]>>value; fprintf(output_fp, " %f", value);}
        fprintf(output_fp, "))\n");
      } catch(YAML::RepresentationException& e) {
      }
    }
  }
  fprintf(output_fp, "\n");

  // limb name
  fprintf(output_fp, "     ;; limbs\n");
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string> link_names = limb.second.first;
    if (link_names.size()>0) fprintf(output_fp, "     (setq %s-root-link %s)\n", limb_name.c_str(), link_names[0].c_str());
    if ( link_names.size() > 0 ) {
      fprintf(output_fp, "     (setq %s (list", limb_name.c_str());
      for (unsigned int i=0;i<link_names.size();i++) fprintf(output_fp, " %s", link_names[i].c_str()); fprintf(output_fp, "))\n");
      fprintf(output_fp, "\n");
    }
  }
  fprintf(output_fp, "\n");

  fprintf(output_fp, "     ;; links\n");
  fprintf(output_fp, "     (setq links (list");
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string> link_names = limb.second.first;
    for (unsigned int i=0;i<link_names.size();i++) fprintf(output_fp, " %s", link_names[i].c_str());
  }
  fprintf(output_fp, "))\n");

  fprintf(output_fp, "     ;; joint-list\n");
  fprintf(output_fp, "     (setq joint-list (list");
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string> joint_names = limb.second.first;
    for (unsigned int i=0;i<joint_names.size();i++) fprintf(output_fp, " %s", joint_names[i].c_str());
  }
  fprintf(output_fp, "))\n");
  fprintf(output_fp, "\n");

#if 0
  // all link name
  fprintf(output_fp, "     ;; links\n");
  fprintf(output_fp, "     (setq links (list ");
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "link", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    fprintf(output_fp, "%s ", thisLink->getName());
  }
  fprintf(output_fp, "))\n");
  fprintf(output_fp, "\n");
  // all joint-list name
  fprintf(output_fp, "     ;; joint-list\n");
  fprintf(output_fp, "     (setq joint-list (list ");
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "joint", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "joint");
    fprintf(output_fp, "%s ", thisLink->getName());
  }
  fprintf(output_fp, "))\n");
  fprintf(output_fp, "\n");
#endif

  // init ending
  fprintf(output_fp, "     ;; init-ending\n");
  fprintf(output_fp, "     (send self :init-ending)\n");
  fprintf(output_fp, "\n");

  // bodies
  fprintf(output_fp, "     ;; overwrite bodies to return draw-things links not (send link :bodies)\n");
  fprintf(output_fp, "     (setq bodies (list ");
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "link", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    fprintf(output_fp, "%s ", thisLink->getName());
  }
  fprintf(output_fp, "))\n");
  fprintf(output_fp, "     self)) ;; :init\n");
  //fprintf(output_fp, "  (:bodies () links)\n");
  fprintf(output_fp, "  )\n\n");

  writeGeometry(output_fp, g_dae->getDatabase());

  return 0;
#if 0
  // node
  int nodeElementCount;
  nodeElementCount = g_dae->getDatabase()->getElementCount(NULL, "node", NULL);
  fprintf(stderr, "Number of Nodes %d\n", nodeElementCount);
  for(int currentNode=0;currentNode<nodeElementCount;currentNode++) {
    domNode *thisNode;
    g_dae->getDatabase()->getElement((daeElement**)&thisNode, currentNode, NULL, "node");
    fprintf(stderr, "%d id:%s, name:%s, sid:%s\n", currentNode, thisNode->getId(), thisNode->getName(), thisNode->getSid());
    int geometryElementCount = thisNode->getInstance_geometry_array().getCount();
    for(int currentGeometry=0;currentGeometry<geometryElementCount;currentGeometry++) {
      domInstance_geometry *thisGeometry = thisNode->getInstance_geometry_array()[currentGeometry];
      fprintf(stderr, " %s\n", thisGeometry->getUrl().id().c_str());
    }
  }

  // get number of links
  int linkElementCount;
  linkElementCount = g_dae->getDatabase()->getElementCount(NULL, "link", NULL);
  fprintf(stderr, "Number of Links %d\n", linkElementCount);

  for(int currentLink=0;currentLink<linkElementCount;currentLink++) {
    // get current link
    domLink *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");

    fprintf(stderr, "link[%2d] %-32s, ", currentLink, thisLink->getName());
    fprintf(stderr, "rotate %d, translate %d, contents %d, full %d\n", 
            thisLink->getRotate_array().getCount(),
            thisLink->getTranslate_array().getCount(),
            thisLink->getContents().getCount(),
            thisLink->getAttachment_full_array().getCount()
            );

    domAxis_constraint_Array jointAxis_array;
    int jointCount;
    if ( thisJoint->getPrismatic_array().getCount() > 0 ) {
      jointAxis_array = thisJoint->getPrismatic_array();
      jointCount = thisJoint->getPrismatic_array().getCount();
    } else if ( thisJoint->getRevolute_array().getCount() > 0 ) {
      jointAxis_array = thisJoint->getRevolute_array();
      jointCount = thisJoint->getRevolute_array().getCount();
    }
    for (int currentJoint=0;currentJoint<jointCount;currentJoint++){
      fprintf(stderr, "axis : %2.0f %2.0f %2.0f, ",
              jointAxis_array[currentJoint]->getAxis()->getValue()[0],
              jointAxis_array[currentJoint]->getAxis()->getValue()[1],
              jointAxis_array[currentJoint]->getAxis()->getValue()[2]);
      fprintf(stderr, "limits :");
      if (jointAxis_array[currentJoint]->getLimits()) {
        fprintf(stderr, " %f %f",
                jointAxis_array[currentJoint]->getLimits()->getMin()->getValue(),
                jointAxis_array[currentJoint]->getLimits()->getMax()->getValue());
      }
    }
    fprintf(stderr, "\n");
  }

  exit(0);
#endif
}

