#include <iostream>
#include "remregs.h"
#include "regdefs.h"
#include "utils.h"
#include "robot.h"

using namespace std;

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
  
  regs.set_reg_b(REG8_MODE, 2);
  ext_key();
  regs.set_reg_b(REG8_MODE, 0);
  
  regs.close();
  return 0;
}
