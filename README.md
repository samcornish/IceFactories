# IceFactories
 Codes used in "Rise and fall of sea ice production in the Arctic Ocean's ice factories", Cornish et al. //
First off, I am sorry that these codes are not tailored for use by people other than me. I wrote them during my PhD, and I don't have time to change them now (I don't work in academia).
But a quick guide to use. I used the server_codes to extract and process data from the original model data files.
Then I used the local_codes to do the analysis. Things got complicated because the files for the first 35 ensemble members were different to members 36 to 40, so there are some scripts in both server_codes and local_codes addressing that issue.
If you choose to use anything here, you'll obviously need to change all the loading and saving information.
The code itself all works and hopefully is vaguely intelligible.
