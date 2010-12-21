#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#include "dae.h"
#include "dom/domCOLLADA.h"

#include <fstream>
#include "yaml-cpp/yaml.h"

#include <boost/foreach.hpp>

extern "C" {
#include <qhull/qhull_a.h>
}


daeDocument *g_document;
DAE* g_dae = NULL;

vector<pair<string, string> > g_all_link_names;

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
  std::vector<coordT> points;

  // get mesh
  domMesh *thisMesh = thisGeometry->getMesh();
  int triangleElementCount = (int)(thisMesh->getTriangles_array().getCount());

  fprintf(fp, "(defclass %s\n", thisGeometry->getId());
  fprintf(fp, "  :super body\n");
  fprintf(fp, "  :slots ())\n");
  fprintf(fp, "(defmethod %s\n", thisGeometry->getId());
  fprintf(fp, "  (:init (&key (name))\n");
  fprintf(fp, "         ;; (replace-object self (make-cube 100 100 100))\n");
  fprintf(fp, "         (replace-object self (send self :qhull-faceset))\n");
  fprintf(fp, "         (if name (send self :name name))\n");
  fprintf(fp, "         self)\n");
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
          float a0,a1,a2;
          a0 = thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3);
          a1 = thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3+1);
          a2 = thisMesh->getSource_array()[0]->getFloat_array()->getValue().get(index*3+2);
          fprintf(fp, "         (gl::glVertex3fv (float-vector %f %f %f))\n",
                  1000*a0,1000*a1,1000*a2);
          // store vertex vector to qhull
          points.push_back(a0);
          points.push_back(a1);
          points.push_back(a2);
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
      fprintf(fp, "     ))\n");

      // do qhull
      char qhull_attr[] = "qhull C-0.001";
      int ret = qh_new_qhull (3, points.size()/3, &points[0], 0, qhull_attr, NULL, NULL);
      if ( ! ret ) {
        fprintf(fp, "  (:qhull-faceset ()\n");
        fprintf(fp, "   ;; qhull %d -> %d faces\n", points.size()/3, qh num_facets);
        fprintf(fp, "   (instance faceset :init :faces (list\n");
        // get faces
        facetT *facet;
        vertexT *vertex, **vertexp;
        FORALLfacets {
          fprintf(fp, "    (instance face :init :vertices (list");
          setT *vertices = qh_facet3vertex(facet); // ccw?
          FOREACHvertex_(vertices) {
            fprintf(fp, " (float-vector %f %f %f)", 1000*vertex->point[0], 1000*vertex->point[1], 1000*vertex->point[2]);
          }
          fprintf(fp, "))\n");
          qh_settempfree(&vertices);
        }
        fprintf(fp, ")))\n");
      }
      qh_freeqhull(!qh_ALL);
      int curlong, totlong;    // memory remaining after qh_memfreeshort
      qh_memfreeshort (&curlong, &totlong);    // free short memory and memory allocator
      if (curlong || totlong) {
        fprintf (stderr, "qhull internal warning (user_eg, #1): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
      }

      fprintf(fp, "  )\n\n");
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
domJoint *findJointFromName(const char *jointName) {
  int jointElementCount;
  jointElementCount = g_dae->getDatabase()->getElementCount(NULL, "joint", NULL);
  for(int currentJoint=0;currentJoint<jointElementCount;currentJoint++) {
    domJoint *thisJoint = NULL;
    // get current geometry
    g_dae->getDatabase()->getElement((daeElement**)&thisJoint, currentJoint, NULL, "joint");
    if ( thisJoint == NULL ) {
      fprintf(stderr, "Counld not found joint %s\n", jointName);
      exit(1);
    }
    if ( thisJoint->getSid() == NULL ) {
      fprintf(stderr, "Counld not found Sid of joint %s\n", jointName);
      exit(1);
    }
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
  for(vector<pair<string, string> >::iterator it=g_all_link_names.begin();it!=g_all_link_names.end();it++){
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
  cerr << "writeJoint " << thisJoint->getName() << ", parent = " << parentLink->getName() << ", child = " << childLink->getName()  << endl;
  fprintf(fp, "     (setq %s\n", thisJoint->getName());
  fprintf(fp, "           (instance %s :init\n",
          (thisJoint->getPrismatic_array().getCount()>0)?"linear-joint":"rotational-joint");
  fprintf(fp, "                     :name :%s\n", findLinkName(thisJoint->getName()));
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

void writeKinematics(FILE *fp, domLink::domAttachment_full_Array thisAttachmentArray) {
  for(unsigned int currentAttachment=0;currentAttachment < thisAttachmentArray.getCount();currentAttachment++) {
    domLinkRef thisLink = thisAttachmentArray[currentAttachment]->getLink();
    writeKinematics(fp, thisLink->getAttachment_full_array());

    cerr << "writeKinematics " << thisLink->getName() << ", childlen = " << thisLink->getAttachment_full_array().getCount() << endl;
    for(unsigned int currentAttachment2=0;currentAttachment2 < (unsigned int)(thisLink->getAttachment_full_array().getCount());currentAttachment2++) {
      writeJoint(fp, thisLink->getAttachment_full_array()[currentAttachment2]->getJoint(),
                 thisLink,
                 thisLink->getAttachment_full_array()[currentAttachment2]->getLink()
                 );
    }
  }
}

void writeTranslate(FILE *fp, const char *indent, const char *name, domTranslate_Array thisArray) {
  int translateCount = thisArray.getCount();
  for(int currentTranslate=0;currentTranslate<translateCount;currentTranslate++){
    domTranslateRef thisTranslate = thisArray[currentTranslate];
    if ( thisTranslate->getSid() ) continue;
    fprintf(fp, "%s(send %s :translate (float-vector %f %f %f))\n",
            indent, name,
            1000*thisTranslate->getValue()[0],
            1000*thisTranslate->getValue()[1],
            1000*thisTranslate->getValue()[2]);
  }
}

void writeRotate(FILE *fp, const char *indent, const char *name, domRotate_Array thisArray) {
  int rotateCount = thisArray.getCount();
  for(int currentRotate=0;currentRotate<rotateCount;currentRotate++){
    domRotateRef thisRotate = thisArray[currentRotate];
    if ( thisRotate->getSid() ) continue;
    fprintf(fp, "%s(send %s :rotate %f (float-vector %f %f %f))\n",
            indent, name,
            (thisRotate->getValue()[3])*M_PI/180.0,
            thisRotate->getValue()[0],
            thisRotate->getValue()[1],
            thisRotate->getValue()[2]);
  }
}

void writeNodes(FILE *fp, domNode_Array thisNodeArray) {
  int nodeArrayCount = thisNodeArray.getCount();
  for(int currentNodeArray=0;currentNodeArray<nodeArrayCount;currentNodeArray++) {
    domNode *thisNode = thisNodeArray[currentNodeArray];
    writeNodes(fp, thisNode->getNode_array());

    if ( strcmp(thisNode->getName(),"visual") == 0 ) continue; //@@@ OK??
    // link
    fprintf(stderr, "link sid:%s name:%s node_array:%d\n",
            thisNode->getSid(), thisNode->getName(),thisNode->getNode_array().getCount() );

    // geometry we assume Node_array()[0] contatins geometry
    if ( thisNode->getNode_array().getCount() > 0 &&
         thisNode->getNode_array()[0]->getInstance_geometry_array().getCount() > 0 ) {
      domNode *thisNode2 = thisNode->getNode_array()[0];
      int geometryCount = thisNode2->getInstance_geometry_array().getCount();
      domInstance_geometry *thisGeometry = thisNode2->getInstance_geometry_array()[0];
      const char * geometryName = (string("b_")+thisGeometry->getUrl().id()).c_str();
      assert(geometryCount == 1);
      fprintf(fp, "     ;; define bodyset-link for %s : %s\n", thisNode->getName(), thisNode->getId());
      fprintf(fp, "     (let (%s)\n", geometryName);
      fprintf(fp, "       (setq %s (instance %s :init))\n",  geometryName, thisGeometry->getUrl().id().c_str());
      // translate
      writeTranslate(fp, "       ", geometryName, thisNode2->getTranslate_array());
      // rotate
      writeRotate(fp, "       ", geometryName, thisNode2->getRotate_array());

      // bodyset link
      fprintf(fp, "       (setq %s\n", thisNode->getName());
      fprintf(fp, "             (instance bodyset-link\n");
      fprintf(fp, "                       :init (make-cascoords)\n");
      fprintf(fp, "                       :bodies (list %s)\n", geometryName);
      fprintf(fp, "                       :name :%s))\n", thisNode->getName());

      // assoc
      for(unsigned int currentNodeArray=0;currentNodeArray<thisNode->getNode_array().getCount();currentNodeArray++) {
        if ( strcmp(thisNode->getNode_array()[currentNodeArray]->getName(),"visual") == 0 ) continue; //@@@ OK??
        fprintf(fp, "       (send %s :assoc %s)\n",
                thisNode->getName(),
                thisNode->getNode_array()[currentNodeArray]->getName());
      }
      fprintf(fp, "       )\n");
    } else if ( thisNode->getNode_array().getCount() > 0 &&
                strcmp(thisNode->getNode_array()[0]->getName(),"visual") != 0 ) {
      cerr << ";; WARNING link without geometry : " << thisNode->getName() << endl;
      fprintf(fp, "     ;; define bodyset-link for %s\n", thisNode->getName());
      fprintf(fp, "     (setq %s (instance bodyset-link :init (make-cascoords) :bodies (list (make-cube 10 10 10)) :name :%s))\n", thisNode->getName(), thisNode->getName());

      // assoc
      for(unsigned int currentNodeArray=0;currentNodeArray<thisNode->getNode_array().getCount();currentNodeArray++) {
        if ( strcmp(thisNode->getNode_array()[currentNodeArray]->getName(),"visual") == 0 ) continue; //@@@ OK??
        fprintf(fp, "     (send %s :assoc %s)\n",
                thisNode->getName(),
                thisNode->getNode_array()[currentNodeArray]->getName());
      }
    } else {
      fprintf(fp, "     ;; define cascaded-coords for %s\n", thisNode->getName());
      fprintf(fp, "     (setq %s (make-cascoords :name :%s))\n", thisNode->getName(), thisNode->getName());

      // assoc
      for(unsigned int currentNodeArray=0;currentNodeArray<thisNode->getNode_array().getCount();currentNodeArray++) {
        if ( strcmp(thisNode->getNode_array()[currentNodeArray]->getName(),"visual") == 0 ) continue; //@@@ OK??
        fprintf(fp, "     (send %s :assoc %s)\n",
                thisNode->getName(),
                thisNode->getNode_array()[currentNodeArray]->getName());
      }
    }

    // translate
    writeTranslate(fp, "     ", thisNode->getName(), thisNode->getTranslate_array());
    // rotate
    writeRotate(fp, "     ", thisNode->getName(), thisNode->getRotate_array());

    fprintf(fp, "\n");
  }
}

int main(int argc, char* argv[]){
  FILE *output_fp;
  char *input_filename, *yaml_filename, *output_filename;
  if (argc!=4) {
    fprintf(stderr, "Usage: %s <input dae filename> <input yaml filename> <output lisp filename>\n",argv[0]);
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
    = {  link_joint_pair("torso", link_joint(torso_link_names, torso_joint_names)),
	 link_joint_pair("larm", link_joint(larm_link_names, larm_joint_names)),
         link_joint_pair("rarm", link_joint(rarm_link_names, rarm_joint_names)),
         link_joint_pair("lleg", link_joint(lleg_link_names, lleg_joint_names)),
         link_joint_pair("rleg", link_joint(rleg_link_names, rleg_joint_names)),
         link_joint_pair("head", link_joint(head_link_names, head_joint_names))
         };

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
          joint_names.push_back(key);
          link_names.push_back(findChildLinkFromJointName(key.c_str())->getName());
          g_all_link_names.push_back(pair<string, string>(key, value));
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
  g_dae->getDatabase()->getElement((daeElement**)&thisVisualscene, 0, NULL, "visual_scene");
  int nodeCount = thisVisualscene->getNode_array().getCount();
  fprintf(stderr, "Number of Nodes %d (= 1)\n", nodeCount); // this shoule be 1
  domNode* thisNode= thisVisualscene->getNode_array()[0];

  fprintf(stderr, "Visual_scene %s\n", thisNode->getName());

  fprintf(output_fp, ";;\n");
  fprintf(output_fp, ";; DO NOT EDIT THIS FILE\n");
  fprintf(output_fp, ";;\n");
  fprintf(output_fp, ";; this file is automatically generated from %s\n", input_filename);
  fprintf(output_fp, ";;\n");
  fprintf(output_fp, ";; %s $ ", get_current_dir_name());for(int i=0;i<argc;i++) fprintf(output_fp, "%s ", argv[i]); fprintf(output_fp, "\n");
  fprintf(output_fp, ";;\n");
  fprintf(output_fp, "\n");
  fprintf(output_fp, "(defun %s () (setq *%s* (instance %s-robot :init)))\n", thisNode->getName(), thisNode->getName(), thisNode->getName());
  fprintf(output_fp, "\n");
  fprintf(output_fp, "(defclass %s-robot\n", thisNode->getName());
  fprintf(output_fp, "  :super robot-model\n");
  fprintf(output_fp, "  :slots (");
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
  fprintf(output_fp, "))\n");

  fprintf(output_fp, "(defmethod %s-robot\n", thisNode->getName());
  fprintf(output_fp, "  (:init\n");
  fprintf(output_fp, "   (&rest args)\n");
  fprintf(output_fp, "   (let ()\n");

  // send super :init
  fprintf(output_fp, "     (send-super* :init :name \"%s\" args)\n", thisNode->getName());
  fprintf(output_fp, "\n");

  // write kinemtaics
  writeNodes(output_fp, thisNode->getNode_array());
  fprintf(output_fp, "     (send self :assoc %s)\n", thisNode->getNode_array()[0]->getName());

  // write joint
  domKinematics_model *thisKinematics;
  g_dae->getDatabase()->getElement((daeElement**)&thisKinematics, 0, NULL, "kinematics_model");
  writeKinematics(output_fp, thisKinematics->getTechnique_common()->getLink_array()[0]->getAttachment_full_array());

  domLinkRef thisLink = thisKinematics->getTechnique_common()->getLink_array()[0];
  for(unsigned int currentAttachment2=0;currentAttachment2 < (unsigned int)(thisLink->getAttachment_full_array().getCount());currentAttachment2++) {
    writeJoint(output_fp, thisLink->getAttachment_full_array()[currentAttachment2]->getJoint(),
               thisLink,
               thisLink->getAttachment_full_array()[currentAttachment2]->getLink()
               );
  }

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
      fprintf(output_fp, "     (send %s :assoc %s-end-coords)\n", link_names.back().c_str(), limb_name.c_str());
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

  domNode *rootNode = thisNode->getNode_array()[0];
  while (rootNode->getNode_array()[0]->getInstance_geometry_array().getCount()==0) {
    rootNode = rootNode->getNode_array()[0];
  }
  fprintf(output_fp, "     (setq links (list %s", rootNode->getName());
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    string limb_name = limb.first;
    vector<string> link_names = limb.second.first;
    for (unsigned int i=0;i<link_names.size();i++) fprintf(output_fp, " %s", link_names[i].c_str());
  }
  fprintf(output_fp, "))\n");

  fprintf(output_fp, "     ;; joint-list\n");
  fprintf(output_fp, "     (setq joint-list (list");
  BOOST_FOREACH(link_joint_pair& limb, limbs) {
    vector<string> joint_names = limb.second.second;
    for (unsigned int i=0;i<joint_names.size();i++) fprintf(output_fp, " %s", joint_names[i].c_str());
  }
  fprintf(output_fp, "))\n");
  fprintf(output_fp, "\n");

  // init ending
  fprintf(output_fp, "     ;; init-ending\n");
  fprintf(output_fp, "     (send self :init-ending)\n");
  fprintf(output_fp, "\n");

  // bodies
  fprintf(output_fp, "     ;; overwrite bodies to return draw-things links not (send link :bodies)\n");

  fprintf(output_fp, "     (setq bodies (flatten (mapcar #'(lambda (b) (if (find-method b :bodies) (send b :bodies))) (list");
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "link", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    fprintf(output_fp, " %s", thisLink->getName());
  }
  fprintf(output_fp, "))))\n\n");

  fprintf(output_fp, "     (send self :reset-pose) ;; :set reset-pose\n\n");
  fprintf(output_fp, "     self)) ;; :init\n\n");

  try {
    const YAML::Node& n = doc["angle-vector"];
    if ( n.size() > 0 ) fprintf(output_fp, "    ;; pre-defined pose methods\n");
    for(YAML::Iterator it=n.begin();it!=n.end();it++) {
      string name; it.first() >> name;
      fprintf(output_fp, "    (:%s () (send self :angle-vector (float-vector", name.c_str());
      const YAML::Node& v = it.second();
      for(unsigned int i=0;i<v.size();i++){
        double d; v[i] >> d;
        fprintf(output_fp, " %f", d);
      }
      fprintf(output_fp, ")))\n");
    }
  } catch(YAML::RepresentationException& e) {
  }

  // all joint and link name
  for(int currentJoint=0;currentJoint<(int)(g_dae->getDatabase()->getElementCount(NULL, "joint", NULL));currentJoint++) {
    domJoint *thisJoint;
    g_dae->getDatabase()->getElement((daeElement**)&thisJoint, currentJoint, NULL, "joint");
    fprintf(output_fp, "    (:%s (&rest args) (forward-message-to %s args))\n", thisJoint->getName(), thisJoint->getName());
  }
  for(int currentLink=0;currentLink<(int)(g_dae->getDatabase()->getElementCount(NULL, "link", NULL));currentLink++) {
    domJoint *thisLink;
    g_dae->getDatabase()->getElement((daeElement**)&thisLink, currentLink, NULL, "link");
    fprintf(output_fp, "    (:%s (&rest args) (forward-message-to %s args))\n", thisLink->getName(), thisLink->getName());
  }

  fprintf(output_fp, "  )\n\n");

  writeGeometry(output_fp, g_dae->getDatabase());

  fprintf(output_fp, "\n\n(provide :%s \"%s\")\n\n", thisNode->getName(), get_current_dir_name());

  return 0;
}

