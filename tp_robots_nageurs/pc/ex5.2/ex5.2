#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include "utils.h"
#include <math.h>

using namespace std;

const uint8_t RADIO_CHANNEL = 201;         ///< robot radio channel
const char* INTERFACE = "COM1";            ///< robot radio interface

int main()
{
  CRemoteRegs regs;
  float freq, ampl;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  
  cout << "Insert frequency: ";
  cin >> freq;
  cout << "Insert amplitude: ";
  cin >> ampl;

  cout << freq << " " << ampl << endl;

  regs.set_reg_b(REG8_MODE, 2);

  regs.set_reg_b(10, ENCODE_PARAM_8(freq,(0),(2)));
  regs.set_reg_b(11, ENCODE_PARAM_8(ampl,(-60),(60)));

  ext_key();

  regs.set_reg_b(REG8_MODE, 0);

  regs.close();
  return 0;
}
