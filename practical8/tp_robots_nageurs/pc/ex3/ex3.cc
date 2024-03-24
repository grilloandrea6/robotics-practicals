#include <iostream>
#include "remregs.h"
#include "robot.h"

using namespace std;

const uint8_t RADIO_CHANNEL = 201;         ///< robot radio channel
const char* INTERFACE = "COM1";            ///< robot radio interface

// Displays the contents of a multibyte register as a list of bytes
void display_multibyte_register(CRemoteRegs& regs, const uint8_t addr)
{
  uint8_t data_buffer[32], len;
  if (regs.get_reg_mb(addr, data_buffer, len) && len == 4) {
    for(uint8_t i = 0; i < 4; i++) {
      cout << (int) (int8_t) data_buffer[i] << " ";
    }
    cout << endl;
  } else {
    cerr << "Data transmission error." << endl;
  }
}

int main()
{
  CRemoteRegs regs;

  if (!init_radio_interface(INTERFACE, RADIO_CHANNEL, regs)) {
    return 1;
  }

  // Reboots the head microcontroller to make sure it is always in the same state
  reboot_head(regs);
  
  while(1) {
    display_multibyte_register(regs, 2);
  }
  
  regs.close();
  return 0;
}
