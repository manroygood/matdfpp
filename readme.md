# New PPLANE and DFIELD

This contains a new front end and many, many fixes to **dfield** and **pplane**, that make them compatible with MATLAB Release 2020b. To distinguish these new working programs, I've renamed them

**matdfield** draws direction fields for (possibly nonautonomous) first-order scalar ordinary differential equations.

**matpplane** draws phase plane portraits for autonomous first-order systems of two ordinary differential equations.

I've also created a simple launcher packaged as a *MATLAB App*, so all a student needs to do is click a button.

There are lots of nice features: 

* drawing level sets
* finding equilibria and linearizing about them.
* finding stable & unstable manifolds

These programs were originally copyrighted by John Polking between 1995 and 2003. The textbook Ordinary Differential Equations Using MATLAB 3rd edition contains a manual for the programs, but they are pretty self-explanatory.

The codes were last compatible with MATLAB 6.5 which was released in 2003. On the MATLAB File Exchange, you can find many submissions that get them running again, but none that attempted to fix all the broken components. 

I have made some attempt to fix all the errors but have not managed to fix everything. Using the Code Compatibility Report and the Code Analyzer Report features in MATLAB, I have examined every warning thrown up by the editor. I have been able to fix almost all of these. Where I could not, I suppressed the warning and left a comment. Most of these are arrays that grow over iterations, and there's not much to do about them. Many came from an overuse of the `eval` command, which I have replaced by more modern coding constructs.

As far as I can tell everything works except for:

* The 'Revert' button 
* The 'Zoom Back' feature in pplane (works in dfield, should be able to fix it pplane).
* "Find a nearly closed orbit" in pplane.
* The message section at the bottom of the display window works great in dfield but is garbled in pplane. 
* Certain forms of complicated ODE cause crashes.

The authors managed to program a large number of useful features but under the hood, the program is a mess. *There is a lot I would like to fix, but the structure of the original function would make modernization nearly impossible.* 

* Each program (pplane & dfield) uses a structure where each callback, instead of being a separate function, is given by a recursive call to the main program (matdfield) or (matpplane), determined by a very long `if-then-else` structure. I have tried breaking out the contents of individual `elseif` blocks into separate functions, in their own m-files. This has generally led to strange and hard-to-trace errors, so I abandoned the effort. This structure makes understanding the scope of variables very challenging.
* The function of each `elseif` block is rarely documented, although I added some documentation in my attempt to understand what each does.
* A lot of code is dedicated to fiddling with window dimensions and similar display issues. The variable names used to do this are not self-explanatory and contain lots of magic numbers.
* A lot of code is dedicated to string manipulation. Either this was written before MATLAB added straightforward string-handling features, or else the authors didn't know about these features. To be fair, the original programs were written before Google existed(!), so discovering MATLAB functionality involved searching in books!
* There are a lot of calls to `eval` and `feval` which leads to hard-to-decipher codes. I have corrected this in the cases where such constructions confused the Code Analyzer.

Nonetheless, there are a lot of simple cleanups that might make this more maintainable, including lots of repeated code that should be moved to separate m-files and made into functions. I will continue to maintain this code and make sure that it remains compatible with future MATLAB releases. If anyone would like to help out, or finds an error they'd like to fix, the code will be available at https://github.com/manroygood. Any other bugs you find, please let me know.

**Acknowledgments**: The current codes are based on revisions by Nathanael Kazmierczak (pplane) and Giampiero Campa (dfield). I gratefully acknowledge their efforts and those of others who posted fixes MATLAB File Exchange. Several other people sent me partially working codes including Joceline Lega, Ross Parker, and Bjorn Sanstede. I got lots of help by posting questions to MATLAB Answers. These were usually answered (within minutes) by Walter Roberson, of the Mathworks. Div Tiwari, also of the Mathworks, provided additional useful suggestions.



