function [A,Bu,Bw,Cz,Dzu,Dzw,Cy,Dyu,Dyw] = model_2fb(rmp,Kg,Kt,Rm,Jm,Km,Beq,Kf1,Kf2,Mf1,Mf2,Mc)

A=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 0 (rmp^2)*Mc*Kf2/(Mc*(rmp^2)*Mf2+Jm*Kg^2*(Mc+Mf2)) -(Mc+Mf2)*(Beq*Rm*(rmp^2)+(Kg^2)*Kt*Km)/(Rm*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) 0 0;
   0 -Kf1/Mf1 Kf2/Mf1 0 0 0;
   0 Kf1/Mf1 -Kf2*(Mc*(rmp^2)*(Mf1+Mf2)+Jm*(Kg^2)*(Mc+Mf1+Mf2))/(Mf1*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) Mc*(Rm*Beq*(rmp^2)+(Kg^2)*Kt*Km)/(Rm*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) 0 0];

Bu=[0 0 0 Kg*Kt*rmp*(Mc+Mf2)/(Rm*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) 0 -Mc*rmp*Kt*Kg/(Rm*((Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))))]';

Bw=[0 0 0 0 -1 0]';

Cz = [0 1 0 0 0 0;0 0 1 0 0 0];

Dzu = [0 0]';

Dzw = [0 0]';

Cy=[1 0 0 0 0 0;
   0 -Kf1/Mf1 Kf2/Mf1 0 0 0;
   0 Kf1/Mf1 -Kf2*(Mc*(rmp^2)*(Mf1+Mf2)+Jm*(Kg^2)*(Mc+Mf1+Mf2))/(Mf1*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) Mc*(Rm*Beq*(rmp^2)+(Kg^2)*Kt*Km)/(Rm*(Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))) 0 0];

Dyu=[0 0 -Mc*rmp*Kt*Kg/(Rm*((Mc*(rmp^2)*Mf2+Jm*(Kg^2)*(Mc+Mf2))))]';

Dyw=[0 -1 0]';



end