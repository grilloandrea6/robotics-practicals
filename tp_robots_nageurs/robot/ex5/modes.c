#include "config.h"
#include "modes.h"
#include "robot.h"
#include "module.h"
#include "registers.h"
#include "hardware.h"

const float FREQ = 1.0;   // Hz
const uint8_t MOTOR_ADDR = 21;

void sine_demo_mode()
{
  uint32_t dt, cycletimer;
  float my_time, delta_t, l;
  int8_t l_rounded;

  cycletimer = getSysTICs();
  my_time = 0;

  init_body_module(MOTOR_ADDR);
  start_pid(MOTOR_ADDR);

  do {
    // Calculates the delta_t in seconds and adds it to the current time
    dt = getElapsedSysTICs(cycletimer);
    cycletimer = getSysTICs();
    delta_t = (float) dt / sysTICSperSEC;
    my_time += delta_t;

    // Calculates the sine wave
    l = 40 * sin(M_TWOPI * FREQ * my_time);
    l_rounded = (int8_t) l;

    bus_set(MOTOR_ADDR, MREG_SETPOINT, DEG_TO_OUTPUT_BODY(l_rounded));
    
    // Make sure there is some delay, so that the timer output is not zero
    pause(ONE_MS);

  } while (reg8_table[REG8_MODE] == IMODE_SINE_DEMO);

  // Back to the "normal" green
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
      case IMODE_SINE_DEMO:
        sine_demo_mode();
        break;
      default:
        reg8_table[REG8_MODE] = IMODE_IDLE;
    }
  }
}