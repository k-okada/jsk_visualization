////////////////////////////////////////////////////////////////////////////////
//
// JSK CR(ColorRange) Capture
//

#include <ros/ros.h>
#include <ros/names.h>

#include <sensor_msgs/Image.h>
#include <sensor_msgs/CameraInfo.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/point_cloud_conversion.h>

#include <image_transport/image_transport.h>
#include <image_transport/subscriber_filter.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

#include <opencv/cv.h>
#include <cv_bridge/CvBridge.h>

#include <tf/transform_listener.h>

#include <dynamic_reconfigure/server.h>
#include "cr_capture/CRCaptureConfig.h"

class CRCaptureNode {
private:
  ros::NodeHandle nh_;
  image_transport::ImageTransport it_;
  image_transport::CameraSubscriber camera_sub_l_, camera_sub_r_, camera_sub_depth_;
  tf::TransformListener tf_;
  ros::Publisher cloud_pub_;
  ros::Publisher cloud2_pub_;
  std::string left_ns_, right_ns_, range_ns_;

  // all subscriber
  message_filters::Subscriber<sensor_msgs::Image> image_conf_sub_;
  message_filters::Subscriber<sensor_msgs::Image> image_intent_sub_;
  message_filters::Subscriber<sensor_msgs::Image> image_depth_sub_;
  message_filters::Subscriber<sensor_msgs::CameraInfo> info_sub_;
  message_filters::TimeSynchronizer<sensor_msgs::Image, sensor_msgs::Image,
				    sensor_msgs::Image, sensor_msgs::CameraInfo> sync_;

  // parameter
  double max_range;
  bool calc_pixelpos;
  bool use_filter;
  int filter_type;
  double trans_pos[3];
  double trans_quat[4];
  bool use_images;
  int intensity_threshold, confidence_threshold;
  bool clear_uncolored_points;
  bool use_smooth;
  bool short_range;
  int smooth_size;
  double smooth_depth,smooth_space;
  double edge1, edge2;
  int dilate_times;

  // dynamic reconfigure
  typedef dynamic_reconfigure::Server<cr_capture::CRCaptureConfig> ReconfigureServer;
  ReconfigureServer reconfigure_server_;

  // buffers
  //sensor_msgs::PointCloud pts_;
  sensor_msgs::CameraInfo info_left_, info_right_, info_depth_;
  IplImage *ipl_left_, *ipl_right_, *ipl_depth_;
  float *map_x, *map_y, *map_z;
  int srheight, srwidth;
  CvMat *cam_matrix, *dist_coeff;
  tf::StampedTransform cam_trans_;

public:
  CRCaptureNode () : nh_("~"), it_(nh_),
		     tf_(ros::Duration(30.0), true), sync_(5),
		     map_x(0), map_y(0), map_z(0) {
    // initialize
    ipl_left_ = new IplImage();
    ipl_right_ = new IplImage();
    ipl_depth_ = new IplImage();

    cam_matrix = cvCreateMat(3, 3, CV_64F);
    dist_coeff = cvCreateMat(1, 5, CV_64F);
    cvSetZero(cam_matrix);
    cvSetZero(dist_coeff);
    cvmSet(cam_matrix, 2, 2, 1.0);

    // Set up dynamic reconfiguration
    ReconfigureServer::CallbackType f = boost::bind(&CRCaptureNode::config_cb, this, _1, _2);
    reconfigure_server_.setCallback(f);

    // parameter
    nh_.param("max_range", max_range, 5.0);
    ROS_INFO("max_range : %f", max_range);

    nh_.param("calculation_pixel", calc_pixelpos, false); // not using
    ROS_INFO("calculation_pixel : %d", calc_pixelpos);

    trans_pos[0] = trans_pos[1] = trans_pos[2] = 0;
    if (nh_.hasParam("translation")) {
      XmlRpc::XmlRpcValue param_val;
      nh_.getParam("translation", param_val);
      if (param_val.getType() == XmlRpc::XmlRpcValue::TypeArray && param_val.size() == 3) {
        trans_pos[0] = param_val[0];
        trans_pos[1] = param_val[1];
        trans_pos[2] = param_val[2];
      }
    }
    ROS_INFO("translation : [%f, %f, %f]", trans_pos[0], trans_pos[1], trans_pos[2]);

    trans_quat[0] = trans_quat[1] = trans_quat[2] = 0;
    trans_quat[3] = 1;
    if (nh_.hasParam("rotation")) {
      XmlRpc::XmlRpcValue param_val;
      nh_.getParam("rotation", param_val);
      if (param_val.getType() == XmlRpc::XmlRpcValue::TypeArray && param_val.size() == 4) {
        trans_quat[0] = param_val[0];
        trans_quat[1] = param_val[1];
        trans_quat[2] = param_val[2];
        trans_quat[3] = param_val[3];
      }
    }
    ROS_INFO("rotation : [%f, %f, %f, %f]", trans_quat[0], trans_quat[1],
             trans_quat[2], trans_quat[3]);

    nh_.param("use_filter", use_filter, true);
    ROS_INFO("use_filter : %d", use_filter);
    edge1 = 40.0; edge2 = 80; dilate_times = 1;

    nh_.param("use_smooth", use_smooth, false);
    ROS_INFO("use_smooth : %d", use_smooth);
    if (use_smooth) {
      nh_.param("smooth_size", smooth_size, 6);
      ROS_INFO("smooth_size : %d", smooth_size);
      nh_.param("smooth_depth", smooth_depth, 0.04);
      ROS_INFO("smooth_depth : %f", smooth_depth);
      smooth_depth = (smooth_depth / max_range) * 0xFFFF;
      nh_.param("smooth_space", smooth_space, 6.0);
      ROS_INFO("smooth_space : %f", smooth_space);
    }

    nh_.param("clear_uncolored_points", clear_uncolored_points, true);
    ROS_INFO("clear_uncolored_points : %d", clear_uncolored_points);

    nh_.param("short_range", short_range, false);
    ROS_INFO("short_range : %d", short_range);

    // ros node setting
    //cloud_pub_ = nh_.advertise<sensor_msgs::PointCloud> ("color_pcloud", 1, msg_connect, msg_disconnect);
    cloud_pub_ = nh_.advertise<sensor_msgs::PointCloud> ("color_pcloud", 1);
    cloud2_pub_ = nh_.advertise<sensor_msgs::PointCloud2> ("color_pcloud2", 1);

    left_ns_ = nh_.resolveName("left");
    right_ns_ = nh_.resolveName("right");
    range_ns_ = nh_.resolveName("range");

    camera_sub_l_ = it_.subscribeCamera(left_ns_ + "/image", 1, &CRCaptureNode::cameraleftCB, this);
    camera_sub_r_ = it_.subscribeCamera(right_ns_ + "/image", 1, &CRCaptureNode::camerarightCB, this);

    // use images for confidence, intensity threshold
    nh_.param("use_images", use_images, false);
    ROS_INFO("use_images : %d", use_images);
    if(use_images) {
      // all subscribe
      nh_.param("intensity_threshold", intensity_threshold, -1);
      if(intensity_threshold >= 0) {
	ROS_INFO("intensity_threshold : %d", intensity_threshold);
      }
      nh_.param("confidence_threshold", confidence_threshold, -1);
      if(confidence_threshold >= 0) {
	ROS_INFO("confidence_threshold : %d", confidence_threshold);
      }

      image_conf_sub_.subscribe(nh_, range_ns_ + "/confidence/image_raw", 1);
      image_intent_sub_.subscribe(nh_, range_ns_ + "/intensity/image_raw", 1);
      image_depth_sub_.subscribe(nh_, range_ns_ + "/distance/image_raw16", 1);
      info_sub_.subscribe(nh_, range_ns_ + "/camera_info", 1);

      sync_.connectInput(image_conf_sub_, image_intent_sub_, image_depth_sub_, info_sub_);
      sync_.registerCallback(boost::bind(&CRCaptureNode::cameraallCB, this, _1, _2, _3, _4));
    } else {
      // not all subscribe
      camera_sub_depth_ = it_.subscribeCamera(range_ns_ + "/image", 1, &CRCaptureNode::camerarangeCB, this);
    }
  }

  void config_cb(cr_capture::CRCaptureConfig &config, uint32_t level) {
    //ROS_INFO("config_");
    clear_uncolored_points = config.clear_uncolored;
    short_range = config.short_range;

    use_filter   = config.use_filter;
    edge1        = config.canny_parameter1;
    edge2        = config.canny_parameter2;
    dilate_times = config.dilate_times;

    use_smooth   = config.use_smooth;
    smooth_size  = config.smooth_size;
    smooth_space = config.smooth_space;
    smooth_depth = config.smooth_depth;
    smooth_depth = (smooth_depth / max_range) * 0xFFFF;

    intensity_threshold = config.intensity_threshold;
    confidence_threshold = config.confidence_threshold;
  }

  void cameraleftCB(const sensor_msgs::ImageConstPtr &img,
                    const sensor_msgs::CameraInfoConstPtr &info) {
    //ROS_INFO("%s", __PRETTY_FUNCTION__);
    if((ipl_left_->width != (int)img->width) ||
       (ipl_left_->height != (int)img->height)) {
      ipl_left_ = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 3);
    }
    sensor_msgs::CvBridge bridge;
    cvResize(bridge.imgMsgToCv(img, "rgb8"), ipl_left_);
    info_left_ = *info;
  }
  void camerarightCB(const sensor_msgs::ImageConstPtr &img,
                    const sensor_msgs::CameraInfoConstPtr &info) {
    //ROS_INFO("%s", __PRETTY_FUNCTION__);
    if((ipl_right_->width != (int)img->width) ||
       (ipl_right_->height != (int)img->height)) {
      ipl_right_ = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 3);
    }
    sensor_msgs::CvBridge bridge;
    cvResize(bridge.imgMsgToCv(img, "rgb8"), ipl_right_);
    info_right_ = *info;
  }

  void cameraallCB (const sensor_msgs::ImageConstPtr &img_c,
		    const sensor_msgs::ImageConstPtr &img_i,
		    const sensor_msgs::ImageConstPtr &img_d,
		    const sensor_msgs::CameraInfoConstPtr &info) {
    if((ipl_depth_->width != (int)img_d->width) ||
       (ipl_depth_->height != (int)img_d->height)) {
      ipl_depth_ = cvCreateImage(cvSize(img_d->width, img_d->height), IPL_DEPTH_16U, 1);
      srwidth = img_d->width;
      srheight = img_d->height;
    }
    sensor_msgs::CvBridge bridge;
    cvResize(bridge.imgMsgToCv(img_d), ipl_depth_); // pass through
    info_depth_ = *info;
    if ( (ipl_right_->width != ipl_left_->width) ||
         (ipl_right_->height != ipl_left_->height) ) {
      //ROS_WARN("invalid image");
      return;
    }

    const unsigned char *conf_img = &(img_c->data[0]);
    const unsigned char *intent_img = &(img_i->data[0]);
    unsigned short *depth_img = (unsigned short*)ipl_depth_->imageData;
    int size = img_c->data.size();
    for(int i=0;i<size;i++) {
      if( (conf_img[i] < confidence_threshold) ||
	  (intent_img[i] < intensity_threshold) ) {
	depth_img[i] = 0;
      }
    }
    //
    calculate_color(img_d, info);
  }
  void camerarangeCB(const sensor_msgs::ImageConstPtr &img,
                     const sensor_msgs::CameraInfoConstPtr &info) {
    //ROS_INFO("%s", __PRETTY_FUNCTION__);
    if((ipl_depth_->width != (int)img->width) ||
       (ipl_depth_->height != (int)img->height)) {
      ipl_depth_ = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_16U, 1);
      srwidth = img->width;
      srheight = img->height;
    }
    sensor_msgs::CvBridge bridge;
    cvResize(bridge.imgMsgToCv(img), ipl_depth_); // pass through
    info_depth_ = *info;
    if ( (ipl_right_->width != ipl_left_->width) ||
         (ipl_right_->height != ipl_left_->height) ) {
      //ROS_WARN("invalid image");
      return;
    }
    calculate_color(img, info);
  }

  void calculate_color(const sensor_msgs::ImageConstPtr &img,
		       const sensor_msgs::CameraInfoConstPtr &info) {

    // smooth birateral filter
    if (use_smooth) {
      cv::Mat in_img16(ipl_depth_);
      cv::Mat in_imgf32(ipl_depth_->height, ipl_depth_->width, CV_32FC1);
      cv::Mat out_imgf32(ipl_depth_->height, ipl_depth_->width, CV_32FC1);
      in_img16.convertTo(in_imgf32, CV_32FC1);
      cv::bilateralFilter(in_imgf32, out_imgf32, smooth_size, smooth_depth, smooth_space, cv::BORDER_REPLICATE);
      out_imgf32.convertTo(in_img16, CV_16UC1);
    }

    // filter outlier
    if (use_filter) {
      cv::Mat in_img16(ipl_depth_);
      cv::Mat in_img(in_img16.rows, in_img16.cols, CV_8UC1);
      cv::Mat out_img(in_img16.rows, in_img16.cols, CV_8UC1);
      in_img16.convertTo(in_img, CV_8UC1, 1.0 / ( 1.0 * 256.0));

      cv::Canny(in_img, out_img, edge1, edge2);
      //cv::dilate(out_img, out_img, cv::Mat());
      if(dilate_times >= 1) {
	cv::dilate(out_img, out_img, cv::Mat(), cv::Point(-1, -1), dilate_times);
      }

      unsigned short *sptr = (unsigned short *)in_img16.data;
      unsigned char *cptr = (unsigned char *)out_img.data;
      for(int i=0;i<in_img16.rows*in_img16.cols;i++) {
	if(*cptr++ > 128) {
	  sptr[i] = 0;
	}
      }
    }

    //check transform
    //tf_.lookupTransform("/sr4000_base", "/left_cam_base", info_depth_.header.stamp, cam_trans_);
    btQuaternion btq(trans_quat[0], trans_quat[1], trans_quat[2], trans_quat[3]);
    btVector3 btp(trans_pos[0], trans_pos[1], trans_pos[2]);
    cam_trans_.setOrigin(btp);
    cam_trans_.setRotation(btq);

    // check info and make map
    makeConvertMap();

    sensor_msgs::PointCloud pts_;
    //convert
    pts_.points.resize(srwidth*srheight);
    convert3DPos(pts_);

    // add color
    if(calc_pixelpos) {
      pts_.channels.resize(4);
      pts_.channels[0].name = "rgb";
      pts_.channels[0].values.resize(srwidth*srheight);
      pts_.channels[1].name = "lu";
      pts_.channels[1].values.resize(srwidth*srheight);
      pts_.channels[2].name = "ru";
      pts_.channels[2].values.resize(srwidth*srheight);
      pts_.channels[3].name = "v";
      pts_.channels[3].values.resize(srwidth*srheight);
    } else {
#if 1 // use rgb
      pts_.channels.resize(1);
      pts_.channels[0].name = "rgb";
      pts_.channels[0].values.resize(srwidth*srheight);
#else
      pts_.channels.resize(3);
      pts_.channels[0].name = "r";
      pts_.channels[0].values.resize(srwidth*srheight);
      pts_.channels[1].name = "g";
      pts_.channels[1].values.resize(srwidth*srheight);
      pts_.channels[2].name = "b";
      pts_.channels[2].values.resize(srwidth*srheight);
#endif
    }
    getColorsOfPointsLRCheck(pts_);

    // advertise
    //pts_.header = info->header;
    if (cloud_pub_.getNumSubscribers() > 0) {
      pts_.header = img->header;
      sensor_msgs::PointCloudPtr ptr = boost::make_shared <sensor_msgs::PointCloud> (pts_);
      cloud_pub_.publish(ptr);
    }
    if (cloud2_pub_.getNumSubscribers() > 0) {
      pts_.header = img->header;
      sensor_msgs::PointCloud2 outbuf;
      if(!sensor_msgs::convertPointCloudToPointCloud2 (pts_, outbuf)) {
	ROS_ERROR ("[cr_capture] Conversion from sensor_msgs::PointCloud2 to sensor_msgs::PointCloud failed!");
	return;
      }
      outbuf.width = srwidth;
      outbuf.height = srheight;
      outbuf.row_step = srwidth * outbuf.point_step;
      sensor_msgs::PointCloud2Ptr ptr = boost::make_shared <sensor_msgs::PointCloud2> (outbuf);
      cloud2_pub_.publish(ptr);
    }
  }

  void makeConvertMap () {
    if( (info_depth_.D[0] != cvmGet(dist_coeff, 0, 0)) ||
        (info_depth_.D[1] != cvmGet(dist_coeff, 0, 1)) ||
        (info_depth_.D[2] != cvmGet(dist_coeff, 0, 2)) ||
        (info_depth_.D[3] != cvmGet(dist_coeff, 0, 3)) ||
        (info_depth_.D[4] != cvmGet(dist_coeff, 0, 4)) ||
        (info_depth_.K[3*0 + 0] != cvmGet(cam_matrix, 0, 0)) ||
        (info_depth_.K[3*0 + 2] != cvmGet(cam_matrix, 0, 2)) ||
        (info_depth_.K[3*1 + 1] != cvmGet(cam_matrix, 1, 1)) ||
        (info_depth_.K[3*1 + 2] != cvmGet(cam_matrix, 1, 2)) ) {
      //
      cvmSet(dist_coeff, 0, 0, info_depth_.D[0]);
      cvmSet(dist_coeff, 0, 1, info_depth_.D[1]);
      cvmSet(dist_coeff, 0, 2, info_depth_.D[2]);
      cvmSet(dist_coeff, 0, 3, info_depth_.D[3]);
      cvmSet(dist_coeff, 0, 4, info_depth_.D[4]);
      //
      cvSetZero(cam_matrix);
      cvmSet(cam_matrix, 2, 2, 1.0);
      cvmSet(cam_matrix, 0, 0, info_depth_.K[3*0 + 0]);
      cvmSet(cam_matrix, 0, 2, info_depth_.K[3*0 + 2]);
      cvmSet(cam_matrix, 1, 1, info_depth_.K[3*1 + 1]);
      cvmSet(cam_matrix, 1, 2, info_depth_.K[3*1 + 2]);
      //
      CvMat *src = cvCreateMat(srheight*srwidth, 1, CV_32FC2);
      CvMat *dst = cvCreateMat(srheight*srwidth, 1, CV_32FC2);
      CvPoint2D32f *ptr = (CvPoint2D32f *)src->data.ptr;
      for(int v=0;v<srheight;v++){
        for(int u=0;u<srwidth;u++){
          ptr->x = u;
          ptr->y = v;
          ptr++;
        }
      }
      if(map_x != 0) delete map_x;
      if(map_y != 0) delete map_y;
      if(map_z != 0) delete map_z;
      map_x = new float[srwidth*srheight];
      map_y = new float[srwidth*srheight];
      map_z = new float[srwidth*srheight];

      cvUndistortPoints(src, dst, cam_matrix, dist_coeff, NULL, NULL);
      ptr = (CvPoint2D32f *)dst->data.ptr;
      for(int i=0;i<srheight*srwidth;i++){
        float xx = ptr->x;
        float yy = ptr->y;
        ptr++;
        double norm = sqrt(xx * xx + yy * yy + 1.0);
        map_x[i] = xx / norm;
        map_y[i] = yy / norm;
        map_z[i] = 1.0 / norm;
      }
      cvReleaseMat(&src);
      cvReleaseMat(&dst);
      ROS_INFO("make conversion map");
    } else {
      //ROS_WARN("do nothing!");
    }
  }

  void convert3DPos(sensor_msgs::PointCloud &pts) {
    int lng=(srwidth*srheight);
    unsigned short *ibuf = (unsigned short*)ipl_depth_->imageData;
    geometry_msgs::Point32 *pt = &(pts.points[0]);

    for(int i=0;i<lng;i++){
      double scl = (max_range *  ibuf[i]) / (double)0xFFFF;

      if(short_range) {
	// magic process from calibration results
	if(scl < 1.0) {
	  scl *= ((-7.071927e+02*pow(scl,3) + 1.825608e+03*pow(scl,2) - 1.571370E+03*scl + 1.454931e+03)/1000.0);
	} else if (scl < 1.1) {
	  scl *= ((1000 + (1.9763 / 100.0) * (1100.0 - 1000.0*scl)) /1000.0);
	}
      }

      if(ibuf[i] >= 0xFFF8) { // saturate
	pt->x = 0.0;
	pt->y = 0.0;
	pt->z = 0.0;
      } else {
	pt->x = (map_x[i] * scl);
	pt->y = (map_y[i] * scl);
	pt->z = (map_z[i] * scl);
      }
      pt++;
    }
  }

  void getColorsOfPointsLRCheck(sensor_msgs::PointCloud &pts) {
    geometry_msgs::Point32 *point_ptr = &(pts.points[0]);
    float fx = info_left_.P[0];
    float cx = info_left_.P[2];
    float fy = info_left_.P[5];
    float cy = info_left_.P[6];
    float tr = info_right_.P[3]; // for ROS projection matrix (unit = m)
    //float tr = (info_right_.P[3])/1000.0; // for jsk projection matrix (unit = mm)

    float *lu_ptr = NULL, *ru_ptr = NULL, *v_ptr = NULL;
    if (calc_pixelpos) {
      lu_ptr = &(pts.channels[1].values[0]);
      ru_ptr = &(pts.channels[2].values[0]);
      v_ptr = &(pts.channels[3].values[0]);
    }
    // ROS_INFO("CHECK/ %f %f %f %f - %f", fx, cx, fy, cy, tr);
    unsigned char *imgl = (unsigned char *)ipl_left_->imageData;
    unsigned char *imgr = (unsigned char *)ipl_right_->imageData;
    int w = ipl_left_->width;
    int h = ipl_left_->height;
    int step = ipl_left_->widthStep;

#define getPixel(img_ptr, pix_x, pix_y, color_r, color_g, color_b) \
    { color_r = img_ptr[step*pix_y + pix_x*3 + 0];                 \
      color_g = img_ptr[step*pix_y + pix_x*3 + 1];                 \
      color_b = img_ptr[step*pix_y + pix_x*3 + 2]; }

    int ypos[srwidth];
    int lxpos[srwidth];
    int rxpos[srwidth];
    int col_x[srwidth];
    int lr_use[srwidth];

    // if(calc_pixelpos) for(int i=0; i < 3*srheight*srwidth; i++) pix[i] = -1;

    for(int y=0;y<srheight;y++) {
      for(int x=0;x<srwidth;x++) {
        int index = y*srwidth + x;
        // convert camera coordinates
        btVector3 pos(point_ptr[index].x, point_ptr[index].y, point_ptr[index].z);
        pos = cam_trans_ * pos;
#if 0
        float posx = point_ptr[index].x;
        float posy = point_ptr[index].y;
        float posz = point_ptr[index].z;
#endif
        float posx = pos[0];
        float posy = pos[1];
        float posz = pos[2];
        if(posz > 0.100) { // filtering near points
          lxpos[x] = (int)(fx/posz * posx + cx); // left cam
          rxpos[x] = (int)((fx*posx + tr)/posz + cx); // right cam
          ypos[x] = (int)(fy/posz * posy + cy);
        } else {
          lxpos[x] = -1;
          rxpos[x] = -1;
          ypos[x] = -1;
        }
        //ROS_INFO("%d\t%d\t%d",lxpos[x],rxpos[x],ypos[x]);
      }
      memset(lr_use, 0, sizeof(int)*srwidth);
      memset(col_x, 0x01000000, sizeof(int)*srwidth);

      int max_lx = -1;
      int min_rx = w;
      for(int x=0;x<srwidth;x++) {
        int lx = lxpos[x];
        int ly = ypos[x];

        int pr = srwidth -x -1;
        int rx = rxpos[pr];
        int ry = ypos[pr];

        if ((w > lx ) && (lx >= 0)
            && (h > ly) && (ly >= 0)) {
          if(lx >= max_lx) {
            max_lx = lx;
          } else {
            lr_use[x] = -1; // use right
          }
        }
        if ((w > rx) && (rx >= 0)
            && (h > ry) && (ry >= 0)) {
          if(rx <= min_rx) {
            min_rx = rx;
          } else {
            lr_use[pr] = 1; // use left
          }
        }
      }
      // finding similar color
      unsigned char lcolr=0, lcolg=0, lcolb=0;
      unsigned char rcolr=0, rcolg=0, rcolb=0;
      for(int x=0;x<srwidth;x++) {
        if(lr_use[x]==0) {
          int lx = lxpos[x];
          int rx = rxpos[x];
          int yy = ypos[x];
          if ((w > lx ) && (lx >= 0)
              && (w > rx) && (rx >= 0)
              && (h > yy) && (yy >= 0)) {
            //imgl->getPixel(lx, yy, &lcolr, &lcolg, &lcolb);
            getPixel(imgl, lx, yy, lcolr, lcolg, lcolb);
            //imgr->getPixel(rx, yy, &rcolr, &rcolg, &rcolb);
            getPixel(imgr, rx, yy, rcolr, rcolg, rcolb);

            double norm = 0.0;
            double norm_r = (double)(lcolr - rcolr);
            double norm_g = (double)(lcolg - rcolg);
            double norm_b = (double)(lcolb - rcolb);
            norm += norm_r * norm_r;
            norm += norm_g * norm_g;
            norm += norm_b * norm_b;
            norm = sqrt(norm);

            if(norm < 50.0) { // magic number for the same color
              col_x[x] = ((0xFF & ((lcolr + rcolr)/2)) << 16) |
                ((0xFF & ((lcolg + rcolg)/2)) << 8) |
                (0xFF & ((lcolb + rcolb)/2));
              //if(pix){
              if(calc_pixelpos) {
#if 0
                int pptr = 3*(x + y*srwidth);
                pix[pptr] = lx;
                pix[pptr+1] = rx;
                pix[pptr+2] = yy;
#endif
                int ptr_pos = (x + y*srwidth);
                lu_ptr[ptr_pos] = lx;
                ru_ptr[ptr_pos] = rx;
                v_ptr[ptr_pos]  = yy;
              }
            } else {
              col_x[x] = 0x2000000;
              // find nearest one in next loop
            }
          } else if ((w > lx ) && (lx >= 0)
                     && (h > yy) && (yy >= 0)) {
            // only left camera is viewing
            //imgl->getPixel(lx, yy, &lcolr, &lcolg, &lcolb);
            getPixel(imgl, lx, yy, lcolr, lcolg, lcolb);
            col_x[x] = ((0xFF & lcolr) << 16) | ((0xFF & lcolg) << 8) | (0xFF & lcolb);
#if 0
            if(pix){
              int pptr = 3*(x + y*srwidth);
              pix[pptr] = lx;
              pix[pptr+2] = yy;
            }
#endif
            if(calc_pixelpos) {
              int ptr_pos = (x + y*srwidth);
              lu_ptr[ptr_pos] = lx;
              v_ptr[ptr_pos]  = yy;
            }
          } else if ((w > rx ) && (rx >= 0)
                     && (h > yy) && (yy >= 0)) {
            // only right camera is viewing
            //imgr->getPixel(rx, yy, &rcolr, &rcolg, &rcolb);
            getPixel(imgr, rx, yy, rcolr, rcolg, rcolb);
            col_x[x] = ((0xFF & rcolr) << 16) | ((0xFF & rcolg) << 8) | (0xFF & rcolb);
#if 0
            if(pix){
              int pptr = 3*(x + y*srwidth);
              pix[pptr+1] = rx;
              pix[pptr+2] = yy;
            }
#endif
            if(calc_pixelpos) {
              int ptr_pos = (x + y*srwidth);
              ru_ptr[ptr_pos] = rx;
              v_ptr[ptr_pos]  = yy;
            }
          } else {
	    // did not find corresponding points in image
            col_x[x] = 0xFF0000;
	    if (clear_uncolored_points) {
	      int pidx = y*srwidth + x;
	      point_ptr[pidx].x = 0.0;
	      point_ptr[pidx].y = 0.0;
	      point_ptr[pidx].z = 0.0;
	    }
          }
        } else if (lr_use[x] > 0) {
          // use left
          int lx = lxpos[x];
          int ly = ypos[x];
          //imgl->getPixel(lx, ly, &lcolr, &lcolg, &lcolb);
          getPixel(imgl, lx, ly, lcolr, lcolg, lcolb);
          col_x[x] = ((0xFF & lcolr) << 16) | ((0xFF & lcolg) << 8) | (0xFF & lcolb);
#if 0
          if(pix){
            int pptr = 3*(x + y*srwidth);
            pix[pptr] = lx;
            pix[pptr+2] = ly;
          }
#endif
          if(calc_pixelpos) {
            int ptr_pos = (x + y*srwidth);
            lu_ptr[ptr_pos] = lx;
            v_ptr[ptr_pos]  = ly;
          }
        } else {
          // use right
          int rx = rxpos[x];
          int ry = ypos[x];
          //imgr->getPixel(rx, ry, &rcolr, &rcolg, &rcolb);
          getPixel(imgr, rx, ry, rcolr, rcolg, rcolb);
          col_x[x] = ((0xFF & rcolr) << 16) | ((0xFF & rcolg) << 8) | (0xFF & rcolb);
#if 0
          if(pix){
            int pptr = 3*(x + y*srwidth);
            pix[pptr+1] = rx;
            pix[pptr+2] = ry;
          }
#endif
          if(calc_pixelpos) {
            int ptr_pos = (x + y*srwidth);
            ru_ptr[ptr_pos] = rx;
            v_ptr[ptr_pos]  = ry;
          }
        }
      }
      // checking color of nearest one
      for(int x=0;x<srwidth;x++) {
        if(col_x[x] & 0x02000000) {
          int n = 0x02000000;
          for(int p=0;p<srwidth;p++) {
            if((x+p >= srwidth) &&
               (x-p < 0)) {
              break;
            } else {
              if(x+p < srwidth) {
                if(!(col_x[x+p] & 0xFF000000)){
                  n = col_x[x+p];
                  break;
                }
              }
              if(x-p >= 0) {
                if(!(col_x[x-p] & 0xFF000000)){
                  n = col_x[x-p];
                  break;
                }
              }
            }
          }
          if(!(n & 0xFF000000)) {
            int lx = lxpos[x];
            int rx = rxpos[x];
            int yy = ypos[x];
            //imgl->getPixel(lx, yy, &lcolr, &lcolg, &lcolb);
            getPixel(imgl, lx, yy, lcolr, lcolg, lcolb);
            //imgr->getPixel(rx, yy, &rcolr, &rcolg, &rcolb);
            getPixel(imgr, rx, yy, rcolr, rcolg, rcolb);
            int clr = (n >> 16) & 0xFF;
            int clg = (n >> 8) & 0xFF;
            int clb = (n >> 0) & 0xFF;
            int dif_l = abs(clr - lcolr) + abs(clg - lcolg) + abs(clb - lcolb);
            int dif_r = abs(clr - rcolr) + abs(clg - rcolg) + abs(clb - rcolb);
            if(dif_l < dif_r) {
              col_x[x] = ((0xFF & lcolr) << 16) | ((0xFF & lcolg) << 8) | (0xFF & lcolb);
#if 0
              if(pix){
                int pptr = 3*(x + y*srwidth);
                pix[pptr] = lx;
                pix[pptr+2] = yy;
              }
#endif
              if(calc_pixelpos) {
                int ptr_pos = (x + y*srwidth);
                lu_ptr[ptr_pos] = lx;
                v_ptr[ptr_pos]  = yy;
              }
            } else {
              col_x[x] = ((0xFF & rcolr) << 16) | ((0xFF & rcolg) << 8) | (0xFF & rcolb);
#if 0
              if(pix){
                int pptr = 3*(x + y*srwidth);
                pix[pptr+1] = rx;
                pix[pptr+2] = yy;
              }
#endif
              if(calc_pixelpos) {
                int ptr_pos = (x + y*srwidth);
                ru_ptr[ptr_pos] = rx;
                v_ptr[ptr_pos]  = yy;
              }
            }
          }
        }
      }
      // setting color
#if 1 // use rgb
      float *colv = &(pts.channels[0].values[y*srwidth]);
      for(int x=0;x<srwidth;x++) {
        colv[x] = *reinterpret_cast<float*>(&(col_x[x]));
      }
#else
      float *colr = &(pts.channels[0].values[y*srwidth]);
      float *colg = &(pts.channels[1].values[y*srwidth]);
      float *colb = &(pts.channels[2].values[y*srwidth]);
      for(int x=0;x<srwidth;x++) {
	int col = col_x[x];
        colr[x] = ((col >> 16) & 0xFF) / 255.0;
	colg[x] = ((col >> 8) & 0xFF) / 255.0;
	colb[x] = (col & 0xFF) / 255.0;
      }
#endif
#if 0
      // setting color
      for(int x=0;x<srwidth;x++,tmp_ptr+=3) {
        int col = col_x[x];
        //printf("col = %X, %d\n", col, tmp_ptr);
        colv[tmp_ptr] = ((float)((col >> 16) & 0xFF))/255.0;
        colv[tmp_ptr+1] = ((float)((col >> 8) & 0xFF))/255.0;
        colv[tmp_ptr+2] = ((float)(col & 0xFF))/255.0;
      }
#endif
    } // y_loop
  }
};

int main(int argc, char **argv)
{
  ros::init(argc, argv, "cr_capture");
  //cv::namedWindow(std::string("window"), );
  CRCaptureNode cap_node;

  ros::spin();
  return 0;
}
