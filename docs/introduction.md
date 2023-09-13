# Notes from Erwin on GTAPPY

## Shocks implemented

1. Different population file replacements
2. GDP change
3. Yield changes?

## Example run file

![](img/2022-12-12-10-31-05.png)

At the top we create the project flow object. We will add tasks to this.

## Example executable call

![](img/2022-12-12-10-33-44.png)

### Harder part to do is figuring out how to have the shocks work

Need to create xSets, xSubsets, Shock.  Could use READFROMFILE command in tablo cmf.

## What are the different kinds of shocks

**Uniform** Shockaoall(AGCOM_SM, reg) = uniform 20;

Another more 

aoall(ACTS, REG) = file select from file pointing to a list of all regions. 0 is no shock.

Everything in the exogenous list can be shocks.

Also can SWAP an initially endogenous with exogenous.

E.g. swap aoreg with GDP

What about changing an elasticity?

those are in gtapparm, so write a new .prm (this is not  a shock but is just replacing an input.)

Notice that there are shocks vs data updates.

The elasticities are being calibrated??? if you change the basedata to basedata2, then a different default2.prm, with no shock, the model will replicate the initial database.

If it's a supply response elasticity (as in PNAS) that WILL affect it (unlike above). 

needs to be percentage change over default POP. In base data. Read this in and process it against Eric's file.

Shock pop(REG) = SELECT FROM FILE filename header 4ltr header.
Swap QGDP aoreg
shock aoall (agcom_sm, reg) = select from yield file.
