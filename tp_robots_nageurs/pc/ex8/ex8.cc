#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include "utils.h"

using namespace std;

const uint8_t   RADIO_CHANNEL    = 126;           ///< robot radio channel
const char*     INTERFACE        = "COM3";        ///< robot radio interface
const char*     TRACKING_PC_NAME = "biorobpc6";   ///< host name of the tracking PC
const uint16_t  TRACKING_PORT    = 10502;         ///< port number of the tracking PC
const uint32_t  TEST_DURATION    = 200;            ///< duration of the movement to be tracked

int main(int argc, char* argv[]) // or char** argv 
{
  CRemoteRegs regs;
  float freq = strtod(argv[1], 0),
        ampl = strtod(argv[2], 0),
        phase= strtod(argv[3], 0),
        steering = strtod(argv[4], 0);

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  regs.set_reg_b(REG8_MODE, 3);

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
  }

  regs.set_reg_b(REG8_MODE, 0);

  regs.close();
  return 0;
}