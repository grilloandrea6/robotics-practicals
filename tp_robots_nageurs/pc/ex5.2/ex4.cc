#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "robot.h"
#include <math.h>

using namespace std;

const uint8_t RADIO_CHANNEL = 201;         ///< robot radio channel
const char* INTERFACE = "COM1";            ///< robot radio interface

#define FREQ (1.5f)
#define AMPL (50.0f)

int main()
{
  CRemoteRegs regs;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  
  regs.set_reg_b(REG8_MODE, 1);

  regs.set_reg_b(10, ENCODE_PARAM_8(FREQ,(0),(2)));
  regs.set_reg_b(11, ENCODE_PARAM_8(AMPL,(-60),(60)));

  regs.close();
  return 0;
}
