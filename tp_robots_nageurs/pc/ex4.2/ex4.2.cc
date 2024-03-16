#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include "utils.h"
#include <math.h>

using namespace std;

#define DT .01

const uint8_t RADIO_CHANNEL = 201;         ///< robot radio channel
const char* INTERFACE = "COM1";            ///< robot radio interface

int main()
{
  CRemoteRegs regs;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  
  regs.set_reg_b(REG8_MODE, 1);

  double actTime = time_d();

  while(!kbhit()) {
    while(time_d() - actTime < DT);
    actTime = time_d();

    regs.set_reg_b(10, ENCODE_PARAM_8(40 * sin(2 * M_PI * actTime),-40,40));
  }
  regs.set_reg_b(REG8_MODE, 0);
  
  regs.close();
  return 0;
}
