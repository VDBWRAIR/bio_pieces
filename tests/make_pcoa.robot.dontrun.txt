*** Settings ***
Library         Process
Library         OperatingSystem
Library         Collections 
Suite Teardown          Terminate All Processes    

*** Variables ***
${fasta} =          tests/testinput/short.aln1.fasta

*** Test Cases ***
TestMakePCAReturnCodeIsZero
    ${process_result} =         Run Process     make_pca        ${fasta}
    Log To Console       ${process_result.stdout}
    Log To Console        ${process_result.stderr}
    Should Be Equal As Integers         ${process_result.rc}        0 
    # Check output
    Should Not Contain         ${process_result.stdout}        Error
