function [Phi] = ObjectiveFunction(MSE, TrB)

lamga= 2;
Phi = MSE./((TrB).^((1+lamga)));



