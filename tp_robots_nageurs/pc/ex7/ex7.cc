#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include "utils.h"
#include <stdint.h>
#include <windows.h>
#include "trkcli.h"
#include <math.h>
#include <ctime.h>

using namespace std;

const uint8_t   RADIO_CHANNEL    = 201;           ///< robot radio channel
const char*     INTERFACE        = "COM1";        ///< robot radio interface
const char*     TRACKING_PC_NAME = "biorobpc11";  ///< host name of the tracking PC
const uint16_t  TRACKING_PORT    = 10502;         ///< port number of the tracking PC
const uint8_t   TEST_DURATION    = 10 * 1000;     ///< duration of the movement to be tracked

void saveToCSV(vector<double> x, vector<double> y, string filename);

int main()
{
  CTrackingClient trk;

  // Connects to the tracking server
  if (!trk.connect(TRACKING_PC_NAME, TRACKING_PORT)) {
    return 1;
  }
  
  CRemoteRegs regs;
  float freq, ampl, phase;
  uint32_t frame_time;
  vector<double> x_hist;
  vector<double> y_hist;
  double x,y;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  
  cout << "Insert frequency: ";
  cin >> freq;
  cout << "Insert amplitude: ";
  cin >> ampl;
  cout << "Insert total phase lag: ";
  cin >> phase;

  cout << "Values set: " << freq << ", " << ampl << ", " << phase << endl;

  regs.set_reg_b(10, ENCODE_PARAM_8(freq,(0),(2)));
  regs.set_reg_b(11, ENCODE_PARAM_8(ampl,(0),(60)));
  regs.set_reg_b(12, ENCODE_PARAM_8(ampl,(0.5),(1.5)));

  regs.set_reg_b(REG8_MODE, 3);

  // start time cycle
  time_t startingTime = time(nullptr);
  time_t currentTime = time(nullptr);

  while(currentTime - startingTime < TEST_DURATION)
  {
    // Gets the current position
    if (!trk.update(frame_time)) {
      cerr << "Could not get position from tracking." << endl;
    } else {

      // Gets the ID of the first spot (the tracking system supports multiple ones)
      int id = trk.get_first_id();
        
      // Reads its coordinates (if (id == -1), then no spot is detected)
      if (id != -1 && trk.get_pos(id, x, y)) {
          x_hist.push_back(x);
          y_hist.push_back(y);
      } else {
          cerr << "Could not get position from tracking." << endl;          
      }
    }
    currentTime = time(nullptr);
  }

  regs.set_reg_b(REG8_MODE, 0);

  string time = ctime(&startingTime);

  saveToCSV(x_hist, y_hist, "trajectory-" + time + ".csv");

  float start_x = x_hist.front(), start_y = y_hist.front(), end_x = x_hist.back(), end_y = y_hist.back();

  float dist = sqrt(pow((end_x - start_x),2) + pow((end_y - start_y),2));

  float vel = dist * 1000 / TEST_DURATION;

  cout << "Total distance traveled:\t" << dist << "m" << endl;
  cout << "Average speed:\t\t\t" << vel << "m/s" << endl;

  regs.close();
  return 0;
}


void saveToCSV(vector<double> x, vector<double> y, string filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open the file: " << filename << endl;
        return;
    }
    file << "x,y" << endl;
     
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << "," << y[i] << "\n";
    }

    file.close();
    cout << "Data saved to: " << filename << endl;
}