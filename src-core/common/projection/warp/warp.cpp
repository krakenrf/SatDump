#include "warp_bkd.h"
#include "logger.h"
#include "core/exception.h"
#include <map>
#include "common/utils.h"
#include "resources.h"
#include "core/opencl.h"
#include <chrono>
#include <cmath>

#include "common/geodetic/geodetic_coordinates.h"

namespace satdump
{
    namespace warp
    {
        double lon_shift(double lon, double shift)
        {
            if (shift == 0)
                return lon;
            lon += shift;
            if (lon > 180)
                lon -= 360;
            if (lon < -180)
                lon += 360;
            return lon;
        }

        void shift_latlon_by_lat(double *lat, double *lon, double shift)
        {
            if (shift == 0)
                return;

            double x = cos(*lat * DEG_TO_RAD) * cos(*lon * DEG_TO_RAD);
            double y = cos(*lat * DEG_TO_RAD) * sin(*lon * DEG_TO_RAD);
            double z = sin(*lat * DEG_TO_RAD);

            double theta = shift * DEG_TO_RAD;

            double x2 = x * cos(theta) + z * sin(theta);
            double y2 = y;
            double z2 = z * cos(theta) - x * sin(theta);

            *lon = atan2(y2, x2) * RAD_TO_DEG;
            double hyp = sqrt(x2 * x2 + y2 * y2);
            *lat = atan2(z2, hyp) * RAD_TO_DEG;
        }

        std::shared_ptr<projection::VizGeorefSpline2D> initTPSTransform(WarpOperation &op)
        {
            return initTPSTransform(op.ground_control_points, op.shift_lon, op.shift_lat);
        }

        std::shared_ptr<projection::VizGeorefSpline2D> initTPSTransform(std::vector<projection::GCP> gcps, int shift_lon, int shift_lat)
        {
            std::shared_ptr<projection::VizGeorefSpline2D> spline_transform = std::make_shared<projection::VizGeorefSpline2D>(2);

            // Attach (non-redundant) points to the transformation.
            std::map<std::pair<double, double>, int> oMapPixelLineToIdx;
            std::map<std::pair<double, double>, int> oMapXYToIdx;
            for (int iGCP = 0; iGCP < (int)gcps.size(); iGCP++)
            {
                double final_lon = lon_shift(gcps[iGCP].lon, shift_lon);
                double final_lat = gcps[iGCP].lat;

                shift_latlon_by_lat(&final_lat, &final_lon, shift_lat);

                const double afPL[2] = {gcps[iGCP].x, gcps[iGCP].y};
                const double afXY[2] = {final_lon, final_lat};

                std::map<std::pair<double, double>, int>::iterator oIter(oMapPixelLineToIdx.find(std::pair<double, double>(afPL[0], afPL[1])));

                if (oIter != oMapPixelLineToIdx.end())
                {
                    if (afXY[0] == gcps[oIter->second].lon && afXY[1] == gcps[oIter->second].lat)
                        continue;
                    else
                    {
                        logger->warn("2 GCPs have the same X,Y!");
                        continue;
                    }
                }
                else
                    oMapPixelLineToIdx[std::pair<double, double>(afPL[0], afPL[1])] = iGCP;

                if (oMapXYToIdx.find(std::pair<double, double>(afXY[0], afXY[1])) != oMapXYToIdx.end())
                {
                    logger->warn("2 GCPs have the same Lat,Lon!");
                    continue;
                }
                else
                    oMapXYToIdx[std::pair<double, double>(afXY[0], afXY[1])] = iGCP;

                if (!spline_transform->add_point(afXY[0], afXY[1], afPL))
                {
                    logger->error("Error generating transformer!");
                    // Handle error appropriately
                }
            }

            logger->info("Solving TPS equations for %d GCPs...", gcps.size());
            auto solve_start = std::chrono::system_clock::now();
            bool solved = spline_transform->solve() != 0;
            if (solved)
                logger->info("Solved! Took %f", (std::chrono::system_clock::now() - solve_start).count() / 1e9);
            else
                logger->error("Failure solving!");

            return spline_transform;
        }

        WarpCropSettings choseCropArea(WarpOperation &op)
        {
            WarpCropSettings cset;
            cset.lat_min = -90;
            cset.lat_max = 90;
            cset.lon_min = -180;
            cset.lon_max = 180;
            cset.y_min = 0;
            cset.y_max = op.output_height;
            cset.x_min = 0;
            cset.x_max = op.output_width;

            std::vector<double> lat_values;
            std::vector<double> lon_values;
            for (projection::GCP &g : op.ground_control_points)
            {
                lat_values.push_back(g.lat);
                lon_values.push_back(g.lon);
            }

            double lat_min = 0;
            double lat_max = 0;
            double lon_min = 0;
            double lon_max = 0;
            lat_min = lat_max = avg_overflowless(lat_values);
            lon_min = lon_max = avg_overflowless(lon_values);

            for (projection::GCP &g : op.ground_control_points)
            {
                if (g.lat > lat_max)
                    lat_max = g.lat;
                if (g.lat < lat_min)
                    lat_min = g.lat;

                if (g.lon > lon_max)
                    lon_max = g.lon;
                if (g.lon < lon_min)
                    lon_min = g.lon;
            }

            // Round to integer degrees
            cset.lat_min = floor(lat_min);
            cset.lon_min = floor(lon_min);
            cset.lat_max = ceil(lat_max);
            cset.lon_max = ceil(lon_max);

            if (op.shift_lat == 90)
                cset.lat_max = 90;
            if (op.shift_lat == -90)
                cset.lat_min = -90;

            // Compute to pixels
            cset.y_max = op.output_height - ((90.0f + cset.lat_min) / 180.0f) * op.output_height;
            cset.y_min = op.output_height - ((90.0f + cset.lat_max) / 180.0f) * op.output_height;
            cset.x_min = (cset.lon_min / 360.0f) * op.output_width + (op.output_width / 2);
            cset.x_max = (cset.lon_max / 360.0f) * op.output_width + (op.output_width / 2);

            // Pixels can offset it a bit - recompute to be 100% accurate
            cset.lat_max = ((op.output_height - cset.y_min) / (double)op.output_height) * 180.0f - 90.0f;
            cset.lat_min = ((op.output_height - cset.y_max) / (double)op.output_height) * 180.0f - 90.0f;
            cset.lon_min = (cset.x_min / (double)op.output_width) * 360.0f - 180.0f;
            cset.lon_max = (cset.x_max / (double)op.output_width) * 360.0f - 180.0f;

            return cset;
        }

        void ImageWarper::warpOnCPU(WarpResult &result)
        {
            // Warp the image on CPU
            auto cpu_start = std::chrono::system_clock::now();
            {
#pragma omp parallel for
                for (int64_t xy_ptr = 0; xy_ptr < (int64_t)result.output_image.width() * (int64_t)result.output_image.height(); xy_ptr++)
                {
                    double xx, yy;
                    double xy[2];
                    int x = (xy_ptr % result.output_image.width());
                    int y = (xy_ptr / result.output_image.width());

                    // Map pixel coordinates to latitude and longitude
                    double lat = -((double)(y + crop_set.y_min) / (double)op.output_height) * 180 + 90;
                    double lon = ((double)(x + crop_set.x_min) / (double)op.output_width) * 360 - 180;

                    // Apply TPS transformation
                    shift_latlon_by_lat(&lat, &lon, op.shift_lat);
                    tps->get_point(lon_shift(lon, op.shift_lon), lat, xy);
                    xx = xy[0];
                    yy = xy[1];

                    // Bounds checking
                    if (xx < 0 || yy < 0)
                        continue;

                    if ((int)xx > (int)op.input_image->width() - 1 || (int)yy > (int)op.input_image->height() - 1)
                        continue;

                    // Pixel value retrieval and assignment
                    if (result.output_image.channels() == 4)
                    {
                        if (op.input_image->channels() == 1)
                            for (int c = 0; c < 3; c++)
                                result.output_image.set(c, y * result.output_image.width() + x, op.input_image->get_pixel_bilinear(0, xx, yy));
                        else if (op.input_image->channels() == 3 || op.input_image->channels() == 4)
                            for (int c = 0; c < 3; c++)
                                result.output_image.set(c, y * result.output_image.width() + x, op.input_image->get_pixel_bilinear(c, xx, yy));

                        if (op.input_image->channels() == 4)
                            result.output_image.set(3, y * result.output_image.width() + x, op.input_image->get_pixel_bilinear(3, xx, yy));
                        else
                            result.output_image.set(3, y * result.output_image.width() + x, 4294967295U);
                    }
                    else
                    {
                        for (int c = 0; c < op.input_image->channels(); c++)
                            result.output_image.set(c, y * result.output_image.width() + x, op.input_image->get_pixel_bilinear(c, xx, yy));
                    }
                }
            }
            auto cpu_time = (std::chrono::system_clock::now() - cpu_start);
            logger->debug("CPU Processing Time %f", cpu_time.count() / 1e9);
        }

#ifdef USE_OPENCL
        void ImageWarper::warpOnGPU_fp32(WarpResult &result)
        {
            // Build GPU Kernel
            cl_program warping_program = opencl::buildCLKernel(resources::getResourcePath("opencl/warp_image_thin_plate_spline_fp32.cl"));

            cl_int err = 0;
            auto &context = satdump::opencl::ocl_context;
            auto &device = satdump::opencl::ocl_device;

            // Now, run the actual OpenCL Kernel
            auto gpu_start = std::chrono::system_clock::now();
            {

                // Images
                cl_mem buffer_map = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(uint32_t) * result.output_image.size(), NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_map! Code " + std::to_string(err));

                cl_mem buffer_img = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * op.input_image->size(), NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_img! Code " + std::to_string(err));

                // TPS Data
                cl_mem buffer_tps_x = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * tps->_nof_points, NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_tps_x! Code " + std::to_string(err));

                cl_mem buffer_tps_y = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * tps->_nof_points, NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_tps_y! Code " + std::to_string(err));

                cl_mem buffer_tps_coefs1 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * tps->_nof_eqs, NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_tps_coefs1! Code " + std::to_string(err));

                cl_mem buffer_tps_coefs2 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * tps->_nof_eqs, NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_tps_coefs2! Code " + std::to_string(err));

                // Image settings
                int img_settings[] = {
                    op.output_width,             // map_img_width
                    op.output_height,            // map_img_height
                    (int)op.input_image->width(),  // img_width
                    (int)op.input_image->height(), // img_height
                    op.input_image->channels(),  // source_channels
                    result.output_image.channels(), // target_channels
                    crop_set.y_min,              // crop_y_min
                    crop_set.y_max,              // crop_y_max
                    crop_set.x_min,              // crop_x_min
                    crop_set.x_max,              // crop_x_max
                    op.shift_lon,                // lon_shiftv
                    op.shift_lat                 // lat_shiftv
                };

                cl_mem buffer_img_settings = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int) * 12, NULL, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create buffer_img_settings! Code " + std::to_string(err));

                // Create an OpenCL queue
                cl_command_queue queue = clCreateCommandQueue(context, device, 0, &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create OpenCL queue! Code " + std::to_string(err));

                // Convert and transfer image data to GPU
                uint16_t* input_image_raw_data = static_cast<uint16_t*>(op.input_image->raw_data());
                size_t input_image_size = op.input_image->size();
                std::vector<uint32_t> input_image_data(input_image_size);
                for (size_t i = 0; i < input_image_size; ++i) {
                    input_image_data[i] = static_cast<uint32_t>(input_image_raw_data[i]);
                }

                err = clEnqueueWriteBuffer(queue, buffer_img, CL_TRUE, 0, sizeof(uint32_t) * input_image_size, input_image_data.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_img! Code " + std::to_string(err));

                // Initialize output image data to zero
                size_t output_image_size = result.output_image.size();
                std::vector<uint32_t> output_image_data(output_image_size, 0);

                err = clEnqueueWriteBuffer(queue, buffer_map, CL_TRUE, 0, sizeof(uint32_t) * output_image_size, output_image_data.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_map! Code " + std::to_string(err));

                // Write TPS data to GPU
                // Convert TPS data to float
                std::vector<float> tps_x(tps->_nof_points);
                for (int i = 0; i < tps->_nof_points; ++i) {
                    tps_x[i] = static_cast<float>(tps->x[i]);
                }
                std::vector<float> tps_y(tps->_nof_points);
                for (int i = 0; i < tps->_nof_points; ++i) {
                    tps_y[i] = static_cast<float>(tps->y[i]);
                }
                std::vector<float> tps_coef1(tps->_nof_eqs);
                for (int i = 0; i < tps->_nof_eqs; ++i) {
                    tps_coef1[i] = static_cast<float>(tps->coef[0][i]);
                }
                std::vector<float> tps_coef2(tps->_nof_eqs);
                for (int i = 0; i < tps->_nof_eqs; ++i) {
                    tps_coef2[i] = static_cast<float>(tps->coef[1][i]);
                }
                float tps_x_mean = static_cast<float>(tps->x_mean);
                float tps_y_mean = static_cast<float>(tps->y_mean);

                err = clEnqueueWriteBuffer(queue, buffer_tps_x, CL_TRUE, 0, sizeof(float) * tps->_nof_points, tps_x.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_tps_x! Code " + std::to_string(err));

                err = clEnqueueWriteBuffer(queue, buffer_tps_y, CL_TRUE, 0, sizeof(float) * tps->_nof_points, tps_y.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_tps_y! Code " + std::to_string(err));

                err = clEnqueueWriteBuffer(queue, buffer_tps_coefs1, CL_TRUE, 0, sizeof(float) * tps->_nof_eqs, tps_coef1.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_tps_coefs1! Code " + std::to_string(err));

                err = clEnqueueWriteBuffer(queue, buffer_tps_coefs2, CL_TRUE, 0, sizeof(float) * tps->_nof_eqs, tps_coef2.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_tps_coefs2! Code " + std::to_string(err));

                // Write image settings
                err = clEnqueueWriteBuffer(queue, buffer_img_settings, CL_TRUE, 0, sizeof(int) * 12, img_settings, 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't write to buffer_img_settings! Code " + std::to_string(err));

                // Initialize the kernel
                cl_kernel warping_kernel = clCreateKernel(warping_program, "warp_image_thin_plate_spline", &err);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't create kernel! Code " + std::to_string(err));

                // Set kernel arguments
                err = clSetKernelArg(warping_kernel, 0, sizeof(cl_mem), &buffer_map);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 0! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 1, sizeof(cl_mem), &buffer_img);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 1! Code " + std::to_string(err));

                // Pass scalar values directly
                err = clSetKernelArg(warping_kernel, 2, sizeof(int), &tps->_nof_points);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 2! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 3, sizeof(cl_mem), &buffer_tps_x);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 3! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 4, sizeof(cl_mem), &buffer_tps_y);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 4! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 5, sizeof(cl_mem), &buffer_tps_coefs1);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 5! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 6, sizeof(cl_mem), &buffer_tps_coefs2);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 6! Code " + std::to_string(err));

                // Pass scalar values directly
                err = clSetKernelArg(warping_kernel, 7, sizeof(float), &tps_x_mean);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 7! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 8, sizeof(float), &tps_y_mean);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 8! Code " + std::to_string(err));

                err = clSetKernelArg(warping_kernel, 9, sizeof(cl_mem), &buffer_img_settings);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 9! Code " + std::to_string(err));
                
                float shift_rad = op.shift_lat * DEG_TO_RAD;
                float cos_shift = cos(shift_rad);
                float sin_shift = sin(shift_rad);

                // Set kernel arguments
                err = clSetKernelArg(warping_kernel, 10, sizeof(float), &cos_shift);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 10! Code " + std::to_string(err));
                err = clSetKernelArg(warping_kernel, 11, sizeof(float), &sin_shift);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't set kernel arg 11! Code " + std::to_string(err));
                

                // Calculate the total number of pixels to process
                //size_t global_work_size = (size_t)(crop_set.x_max - crop_set.x_min) * (size_t)(crop_set.y_max - crop_set.y_min);

                // Enqueue the kernel
                //err = clEnqueueNDRangeKernel(queue, warping_kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
                //if (err != CL_SUCCESS)
                //    throw satdump_exception("Couldn't enqueue kernel! Code " + std::to_string(err));
            
            
                // Get proper workload size
                size_t size_wg = 0;
                size_t compute_units = 0;
                clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &size_wg, NULL);
                clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(size_t), &compute_units, NULL);

                logger->debug("Workgroup size %d", size_wg * compute_units);

                // Run the kernel!
                size_t total_wg_size = int(size_wg) * int(compute_units);            
            
                // Calculate total number of pixels to process
                size_t crop_width = (size_t)(crop_set.x_max - crop_set.x_min);
                size_t crop_height = (size_t)(crop_set.y_max - crop_set.y_min);
                size_t total_pixels = crop_width * crop_height;


                size_t local_work_size = total_wg_size; // Experiment with different values
                size_t global_work_size = ((total_pixels + local_work_size - 1) / local_work_size) * local_work_size;
                
                logger->debug("Global Work Size %d", global_work_size);

                // Enqueue the kernel
                err = clEnqueueNDRangeKernel(queue, warping_kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't enqueue kernel! Code " + std::to_string(err));
                //CHECK_CL_ERROR(err, "Couldn't enqueue kernel");



                // Wait for kernel execution to finish
                clFinish(queue);

                // Read the result back to host
                err = clEnqueueReadBuffer(queue, buffer_map, CL_TRUE, 0, sizeof(uint32_t) * output_image_size, output_image_data.data(), 0, NULL, NULL);
                if (err != CL_SUCCESS)
                    throw satdump_exception("Couldn't read buffer_map! Code " + std::to_string(err));

                // Convert output data from uint32_t to uint16_t if necessary
                uint16_t* output_image_raw_data = static_cast<uint16_t*>(result.output_image.raw_data());
                for (size_t i = 0; i < output_image_size; ++i) {
                    output_image_raw_data[i] = static_cast<uint16_t>(output_image_data[i]);
                }

                // Release resources
                clReleaseMemObject(buffer_img);
                clReleaseMemObject(buffer_map);
                clReleaseMemObject(buffer_tps_x);
                clReleaseMemObject(buffer_tps_y);
                clReleaseMemObject(buffer_tps_coefs1);
                clReleaseMemObject(buffer_tps_coefs2);
                clReleaseMemObject(buffer_img_settings);
                clReleaseKernel(warping_kernel);
                //clReleaseProgram(warping_program);
                clReleaseCommandQueue(queue);
            }
            auto gpu_time = (std::chrono::system_clock::now() - gpu_start);
            logger->debug("GPU Processing Time %f", gpu_time.count() / 1e9);
        }

#endif

        void ImageWarper::update(bool skip_tps)
        {
            if (!skip_tps)
                tps = initTPSTransform(op);
            crop_set = choseCropArea(op);
        }

        WarpResult ImageWarper::warp(bool force_double)
        {
            WarpResult result;

            // Prepare the output image with 32 bits per channel
            result.output_image = image::Image(32,
                                               crop_set.x_max - crop_set.x_min, crop_set.y_max - crop_set.y_min,
                                               op.output_rgba ? 4 : op.input_image->channels());
            result.top_left = {0, 0, (double)crop_set.lon_min, (double)crop_set.lat_max};
            result.top_right = {(double)result.output_image.width() - 1, 0, (double)crop_set.lon_max, (double)crop_set.lat_max};
            result.bottom_left = {0, (double)result.output_image.height() - 1, (double)crop_set.lon_min, (double)crop_set.lat_min};
            result.bottom_right = {(double)result.output_image.width() - 1, (double)result.output_image.height() - 1, (double)crop_set.lon_max, (double)crop_set.lat_min};

#ifdef USE_OPENCL
            if (satdump::opencl::useCL())
            {
                try
                {
                    logger->debug("Using GPU! Double precision requested %d", (int)force_double);
                    satdump::opencl::setupOCLContext();
                    warpOnGPU_fp32(result);
                    return result;
                }
                catch (std::runtime_error &e)
                {
                    logger->error("Error warping on GPU: %s", e.what());
                }
            }
#endif

            logger->debug("Using CPU!");
            warpOnCPU(result);

            return result;
        }
    }
}
