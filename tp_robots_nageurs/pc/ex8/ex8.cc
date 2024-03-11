#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include "utils.h"
#include <stdint.h>
#include <windows.h>
#include "trkcli.h"
#include <math.h>
#include <ctime>
#include <string>
#include <string.h>
#include <vector>
#include <fstream>

using namespace std;

const uint8_t   RADIO_CHANNEL    = 126;           ///< robot radio channel
const char*     INTERFACE        = "COM3";        ///< robot radio interface
const char*     TRACKING_PC_NAME = "biorobpc6";   ///< host name of the tracking PC
const uint16_t  TRACKING_PORT    = 10502;         ///< port number of the tracking PC
const uint32_t  TEST_DURATION    = 200;            ///< duration of the movement to be tracked

void saveToCSV(vector<double> x, vector<double> y, vector<double> t, string filename);

int main(int argc, char* argv[]) // or char** argv 
{

  
  CRemoteRegs regs;
  float freq = strtod(argv[1], 0),
        ampl = strtod(argv[2], 0),
        phase= strtod(argv[3], 0),
        steering = strtod(argv[4], 0);
  int do_pid = atoi(argv[5]);

  vector<double> x_hist;
  vector<double> y_hist;
  vector<double> t_hist;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  if(do_pid)
  {
    // Reboots the head microcontroller to make sure it is always in the same state
    reboot_head(regs);

    regs.set_reg_b(REG8_MODE, 3);
    cin >> do_pid;
  } 
  /*  
  cout << "Insert frequency: ";
  cin >> freq;
  cout << "Insert amplitude: ";
  cin >> ampl;
  cout << "Insert total phase lag: ";
  cin >> phase;
  */
  cout << "Values set: " << freq << ", " << ampl << ", " << phase << endl;

  regs.set_reg_b(10, ENCODE_PARAM_8(freq,(0),(2)));
  regs.set_reg_b(11, ENCODE_PARAM_8(ampl,(0),(60)));
  regs.set_reg_b(12, ENCODE_PARAM_8(phase,(0.5),(1.5)));
  regs.set_reg_b(13, ENCODE_PARAM_8(steering,(-30),(+30)));

  uint32_t rgb = (255 << 16) | (255 << 8) | 255;
  regs.set_reg_dw(0, rgb);


  // start time cycle
  double startingTime = time_d();
  double currentTime  = time_d();

  while(currentTime - startingTime < TEST_DURATION)
  {
    currentTime = time_d();

    uint16_t inp = ext_key();
    char a = inp & 0xFF;
    if(a == 'a') steering = 20;
    else if(a == 'w') steering = 0;
    else if(a == 'd') steering = -20;
    cout << a << endl;

    regs.set_reg_b(13, ENCODE_PARAM_8(steering,(-30),(+30)));

    /*
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
          t_hist.push_back(currentTime);
          cout.precision(8);
          cout << "x: " << x << " y: " << y << " t: " << fixed << currentTime << endl;

      } else {
          cerr << "Could not get position from tracking." << endl;          
      }
    }*/
  }

  regs.set_reg_b(REG8_MODE, 0);
/*
  string my_time = std::to_string(time(nullptr));

  saveToCSV(x_hist, y_hist, t_hist, string("trajectory-") + my_time + "_" + std::to_string(phase) + "_" + std::to_string(freq) + "_" + std::to_string(ampl) + "_" + std::to_string(steering) + ".csv");

  double start_x = x_hist.front(), start_y = y_hist.front(), end_x = x_hist.back(), end_y = y_hist.back();

  double dist = sqrt(pow((end_x - start_x),2) + pow((end_y - start_y),2));

  double vel = dist / TEST_DURATION;

  cout << "Total distance traveled:\t" << dist << "m" << endl;
  cout << "Average speed:\t\t\t" << vel << "m/s" << endl;
*/
  regs.close();
  return 0;
}


void saveToCSV(vector<double> x, vector<double> y, vector<double> t, string filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open the file: " << filename << endl;
        return;
    }
    file << "x,y,t" << endl;

    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << "," << y[i] << "," << t[i] << endl;
    }

    file.close();
    cout << "Data saved to: " << filename << endl;
}