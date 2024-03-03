#include "config.h"
#include "modes.h"
#include "robot.h"
#include "module.h"
#include "registers.h"
#include "hardware.h"

const uint8_t MOTOR_ADDR = 21;

double setpoint = 0;

/* Register callback function, handles some new registers on the radio.
 * All these registers are of course completely useless, but it demonstrates how
 * to implement a register callback function, and what it can do.
 */
static int8_t register_handler(uint8_t operation, uint8_t address, RadioData* radio_data)
{
  uint8_t i;
  
  switch (operation)
  {
    case ROP_WRITE_8:
      if (address == 10) {
        setpoint = DECODE_PARAM_8(radio_data->byte,(-50),(50));
        return TRUE;
      }
  }
  return FALSE;
}



void motor_demo_mode()
{
  // Registers the register handler callback function
  radio_add_reg_callback(register_handler);

  init_body_module(MOTOR_ADDR);
  start_pid(MOTOR_ADDR);
  set_color(4);
  while (reg8_table[REG8_MODE] == IMODE_MOTOR_DEMO) {
    bus_set(MOTOR_ADDR, MREG_SETPOINT, DEG_TO_OUTPUT_BODY(setpoint));
    pause(TEN_MS);
  }
  bus_set(MOTOR_ADDR, MREG_SETPOINT, DEG_TO_OUTPUT_BODY(0.0));
  pause(ONE_SEC);
  bus_set(MOTOR_ADDR, MREG_MODE, MODE_IDLE);
  set_color(2);
}

void main_mode_loop()
{
  reg8_table[REG8_MODE] = IMODE_IDLE;

  while (1)
  {
    switch(reg8_table[REG8_MODE])
    {
      case IMODE_IDLE:
        break;
      case IMODE_MOTOR_DEMO:
        motor_demo_mode();
        break;
      default:
        reg8_table[REG8_MODE] = IMODE_IDLE;
    }
  }
}
