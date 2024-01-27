#include<AccelStepper.h>
#define step1 2
#define step2 3
#define dir1 5
#define dir2 6

AccelStepper M1(1, step1, dir1);
AccelStepper M2(1, step2, dir2);
double R = 200/(TWO_PI)*16;
double Theta1 = 0;
double Theta2 = 0;
double Target1 = 0;
double Target2 = 0;

void setup() {
  Serial.begin(115200);
  pinMode(step1, OUTPUT);
  pinMode(step2, OUTPUT);
  pinMode(dir1, OUTPUT);
  pinMode(dir2, OUTPUT);
  pinMode(19, OUTPUT);
  pinMode(20, OUTPUT);
  pinMode(21, OUTPUT);
  
  pinMode(16, OUTPUT);
  pinMode(17, OUTPUT);
  pinMode(18, OUTPUT);

  M1.setAcceleration(10000000);
  M2.setAcceleration(10000000);

  M1.setMaxSpeed(250*16);
  M2.setMaxSpeed(250*16);
  M1.moveTo(Theta1);
  M2.moveTo(Theta2);
  digitalWrite(19, 0);
  digitalWrite(20, 1);
  digitalWrite(21, 0);
  digitalWrite(16, 0);
  digitalWrite(17, 1);
  digitalWrite(18, 0);
}

void loop() {
  if (Serial.available())
  {
    String in = Serial.readStringUntil('\n');
    if (in == "on")
    {
      digitalWrite(21, 1);
      digitalWrite(18, 1);
    }
    else if (in == "off")
    {
      digitalWrite(21, 0);
      digitalWrite(18, 0);
    }
    else
    {
      Theta1 = getValue(in, ',', 0).toDouble()*R;
      Theta2 = getValue(in, ',', 1).toDouble()*R;
      //double sec = getValue(in, ',', 2).toDouble();
      M1.moveTo((Theta1));
      M2.moveTo((Theta2));
      /*M1.setSpeed((Theta1-Target1)*10/sec*R);
      Target1 = Theta1;
      M2.setSpeed((Theta2-Target2)*10/sec*R);
      Target1 = Theta2;*/
    }
  }
  M1.run();
  M2.run();
}
String getValue(String data, char separator, int index)
{
  int found = 0;
  int strIndex[] = {0, -1};
  int maxIndex = data.length() - 1;

  for (int i = 0; i <= maxIndex && found <= index; i++) {
    if (data.charAt(i) == separator || i == maxIndex) {
      found++;
      strIndex[0] = strIndex[1] + 1;
      strIndex[1] = (i == maxIndex) ? i + 1 : i;
    }
  }

  return found > index ? data.substring(strIndex[0], strIndex[1]) : "";
}
