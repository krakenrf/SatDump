#define DEG_TO_RAD (3.14159265359f / 180.0f)
#define RAD_TO_DEG (180.0f / 3.14159265359f)

inline float lon_shift(float lon, float shift) {
    lon += shift;
    lon -= 360.0f * (lon > 180.0f);
    lon += 360.0f * (lon < -180.0f);
    return lon;
}

inline void shift_latlon_by_lat(float *lat, float *lon, float cos_shift, float sin_shift) {
    if (cos_shift == 1.0f && sin_shift == 0.0f)
        return;

    float lat_rad = *lat * DEG_TO_RAD;
    float lon_rad = *lon * DEG_TO_RAD;

    float cos_lat = native_cos(lat_rad);
    float sin_lat = native_sin(lat_rad);
    float cos_lon = native_cos(lon_rad);
    float sin_lon = native_sin(lon_rad);

    float x = cos_lat * cos_lon;
    float y = cos_lat * sin_lon;
    float z = sin_lat;

    float x2 = x * cos_shift + z * sin_shift;
    float y2 = y;
    float z2 = z * cos_shift - x * sin_shift;

    *lon = atan2(y2, x2) * RAD_TO_DEG;
    float hyp = native_sqrt(x2 * x2 + y2 * y2);
    *lat = atan2(z2, hyp) * RAD_TO_DEG;
}


inline float SQ(const float x) { return x * x; }

inline void VizGeorefSpline2DBase_func4(float *res, const float *pxy,
                                        global const float *xr,
                                        global const float *yr) {
  float dist0 = SQ(xr[0] - pxy[0]) + SQ(yr[0] - pxy[1]);
  res[0] = dist0 != 0.0f ? dist0 * native_log(dist0) : 0.0f;
  float dist1 = SQ(xr[1] - pxy[0]) + SQ(yr[1] - pxy[1]);
  res[1] = dist1 != 0.0f ? dist1 * native_log(dist1) : 0.0f;
  float dist2 = SQ(xr[2] - pxy[0]) + SQ(yr[2] - pxy[1]);
  res[2] = dist2 != 0.0f ? dist2 * native_log(dist2) : 0.0f;
  float dist3 = SQ(xr[3] - pxy[0]) + SQ(yr[3] - pxy[1]);
  res[3] = dist3 != 0.0f ? dist3 * native_log(dist3) : 0.0f;
}

inline float VizGeorefSpline2DBase_func(const float x1, const float y1,
                                        const float x2, const float y2) {
  const float dist = SQ(x2 - x1) + SQ(y2 - y1);
  return dist != 0.0f ? dist * native_log(dist) : 0.0f;
}

inline unsigned int get_pixel_bilinear(global const unsigned int *img, const int width,
                                       const int height, const int channels, int cc,
                                       float rx, float ry) {
  int x = (int)rx;
  int y = (int)ry;

  float x_diff = rx - x;
  float y_diff = ry - y;

  size_t index = (y * width + x);
  size_t max_index = width * height;

  unsigned int a = 0, b = 0, c = 0, d = 0;
  float a_a = 1.0f, b_a = 1.0f, c_a = 1.0f, d_a = 1.0f;

  a = img[cc * width * height + index];
  if (channels == 4 && cc != 3)
    a_a = img[3 * width * height + index] / 4294967295.0f; // Max value of uint32_t

  if (index + 1 < max_index) {
    b = img[cc * width * height + index + 1];
    if (channels == 4 && cc != 3) {
      b_a = img[3 * width * height + index + 1] / 4294967295.0f;
      b = (float)b * b_a;
    }
  } else
    return a;

  if (index + width < max_index) {
    c = img[cc * width * height + index + width];
    if (channels == 4 && cc != 3) {
      c_a = img[3 * width * height + index + width] / 4294967295.0f;
      c = (float)c * c_a;
    }
  } else
    return a;

  if (index + width + 1 < max_index) {
    d = img[cc * width * height + index + width + 1];
    if (channels == 4 && cc != 3) {
      d_a = img[3 * width * height + index + width + 1] / 4294967295.0f;
      d = (float)d * d_a;
    }
  } else
    return a;

  if (x == width - 1)
    return a;
  if (y == height - 1)
    return a;

  a = (float)a * a_a;
  float val_f = (float)a * (1.0f - x_diff) * (1.0f - y_diff) +
                (float)b * x_diff * (1.0f - y_diff) +
                (float)c * y_diff * (1.0f - x_diff) +
                (float)d * x_diff * y_diff;
  if (channels == 4 && cc != 3) {
    float denom = a_a * (1.0f - x_diff) * (1.0f - y_diff) +
                  b_a * x_diff * (1.0f - y_diff) +
                  c_a * y_diff * (1.0f - x_diff) +
                  d_a * x_diff * y_diff;
    val_f = val_f / denom;
  }

  unsigned int val = (unsigned int)(val_f + 0.5f);
  if (val > 4294967295U)
    val = 4294967295U;
  return val;
}

void kernel warp_image_thin_plate_spline(
    global unsigned int *map_image,
    global const unsigned int *img,
    int tps_no_points,
    global const float *tps_x,
    global const float *tps_y,
    global const float *tps_coef_1,
    global const float *tps_coef_2,
    float tps_xmean,
    float tps_ymean,
    global const int *img_settings,
    float cos_shift,
    float sin_shift
) {

  int id = get_global_id(0);

  int map_img_width = img_settings[0];
  int map_img_height = img_settings[1];
  int crop_x_min = img_settings[8];
  int crop_y_min = img_settings[6];
  int crop_width = img_settings[9] - img_settings[8];
  int crop_height = img_settings[7] - img_settings[6];
  int img_width = img_settings[2];
  int img_height = img_settings[3];
  int source_channels = img_settings[4];
  int target_channels = img_settings[5];
  int lon_shiftv = img_settings[10];
  int lat_shiftv = img_settings[11];

  size_t n = crop_width * crop_height;

  if (id >= n)
      return;

  int x = id % crop_width;
  int y = id / crop_width;

  // Adjust x and y based on crop offsets
  x += crop_x_min;
  y += crop_y_min;

  // Map pixel coordinates to latitude and longitude
  float lat = -((float)y / (float)map_img_height) * 180.0f + 90.0f;
  float lon = ((float)x / (float)map_img_width) * 360.0f - 180.0f;

  // Apply TPS transformation
  //shift_latlon_by_lat(&lat, &lon, (float)lat_shiftv);
  //const float Pxy[2] = {lon_shift(lon, (float)lon_shiftv) - tps_xmean, lat - tps_ymean};

  shift_latlon_by_lat(&lat, &lon, cos_shift, sin_shift);
  float lon_shifted = lon_shift(lon, (float)lon_shiftv);
  const float Pxy[2] = {lon_shifted - tps_xmean, lat - tps_ymean};



  float vars[2];
  vars[0] = tps_coef_1[0] + tps_coef_1[1] * Pxy[0] + tps_coef_1[2] * Pxy[1];
  vars[1] = tps_coef_2[0] + tps_coef_2[1] * Pxy[0] + tps_coef_2[2] * Pxy[1];

  int r = 0;
  for (; r < (tps_no_points & (~3)); r += 4) {
    float dfTmp[4] = {};
    VizGeorefSpline2DBase_func4(dfTmp, Pxy, &tps_x[r], &tps_y[r]);
    vars[0] += tps_coef_1[r + 3] * dfTmp[0] + tps_coef_1[r + 3 + 1] * dfTmp[1] +
               tps_coef_1[r + 3 + 2] * dfTmp[2] + tps_coef_1[r + 3 + 3] * dfTmp[3];
    vars[1] += tps_coef_2[r + 3] * dfTmp[0] + tps_coef_2[r + 3 + 1] * dfTmp[1] +
               tps_coef_2[r + 3 + 2] * dfTmp[2] + tps_coef_2[r + 3 + 3] * dfTmp[3];
  }
  for (; r < tps_no_points; r++) {
    const float tmp = VizGeorefSpline2DBase_func(Pxy[0], Pxy[1], tps_x[r], tps_y[r]);
    vars[0] += tps_coef_1[r + 3] * tmp;
    vars[1] += tps_coef_2[r + 3] * tmp;
  }

  float xx = vars[0];
  float yy = vars[1];

  // Bounds checking
  if (xx < 0.0f || yy < 0.0f || xx > (float)(img_width - 1) || yy > (float)(img_height - 1))
    return;

  // Compute position in output image
  int output_x = x - crop_x_min;
  int output_y = y - crop_y_min;

  // Pixel value retrieval and assignment
  if (target_channels == 4) {
    if (source_channels == 1)
      for (int c = 0; c < 3; c++)
        map_image[(crop_width * crop_height) * c + output_y * crop_width + output_x] =
            get_pixel_bilinear(img, img_width, img_height, source_channels, 0, xx, yy);
    else if (source_channels == 3 || source_channels == 4)
      for (int c = 0; c < 3; c++)
        map_image[(crop_width * crop_height) * c + output_y * crop_width + output_x] =
            get_pixel_bilinear(img, img_width, img_height, source_channels, c, xx, yy);
    if (source_channels == 4)
      map_image[(crop_width * crop_height) * 3 + output_y * crop_width + output_x] =
          get_pixel_bilinear(img, img_width, img_height, source_channels, 3, xx, yy);
    else
      map_image[(crop_width * crop_height) * 3 + output_y * crop_width + output_x] = 4294967295U; // Max uint32_t
  } else {
    for (int c = 0; c < source_channels; c++)
      map_image[(crop_width * crop_height) * c + output_y * crop_width + output_x] =
          get_pixel_bilinear(img, img_width, img_height, source_channels, c, xx, yy);
  }
}
