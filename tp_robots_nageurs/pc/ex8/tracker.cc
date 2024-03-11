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
const uint32_t  TEST_DURATION    = 40;            ///< duration of the movement to be tracked

void saveToCSV(vector<double> x, vector<double> y, vector<double> t, string filename);

int main(int argc, char* argv[]) // or char** argv 
{
  CTrackingClient trk;
cout << " eccheppalleeeee" << endl;
  // Connects to the tracking server
  if (!trk.connect(TRACKING_PC_NAME, TRACKING_PORT)) {
    return 1;
  }
  
  uint32_t frame_time;
  vector<double> x_hist;
  vector<double> y_hist;
  vector<double> t_hist;
  double x,y;


  // start time cycle
  double startingTime = time_d();
  double currentTime  = time_d();
    cout << " eccheppalleeeee" << endl;
  while(currentTime - startingTime < TEST_DURATION)
  {
    cout << "-";
    currentTime = time_d();

    
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
    }
  }

  string my_time = std::to_string(time(nullptr));

  saveToCSV(x_hist, y_hist, t_hist, string("trajectory-") + my_time + "keyboard.csv");

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